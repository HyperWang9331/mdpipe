#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把补好链的受体收拾成能直接进 pdb2gmx 的样子。

补完 loop 还不算完。晶体结构里常见三类"不是标准氨基酸"的东西，pdb2gmx 一个都
不认，全都报 Residue 'XXX' not found in residue topology database：

  1. 修饰氨基酸：TPO（磷酸苏氨酸）、SEP、PTR、CSO、MSE……方法是降级成母体残基
     THR/SER/TYR/CYS/MET，并把多出来的原子（磷酸基的 P/O1P/O2P/O3P、乙酰基的
     C2A/C1A/OAC）删掉。**磷酸化状态会因此丢失**，要保留就得单独参数化。
  2. 共晶小分子被写进同一条链里（BCL2 的 N2L、STAT3 的 KQV 就是），没有主链
     N/CA/C，删掉不影响多肽链的连续性。做受体 MD 本来也该把它们去掉。
  3. 水、离子：默认删除（--keep-water / --keep-ligands 可以留）。

这个工具只动"残基名和多余原子"，不动坐标、不动主链，所以补链结果不会被改掉。

用法：
    python scripts/md_ready.py rec_modeled.pdb --out rec_mdready.pdb
    python scripts/md_ready.py rec.pdb --out rec_mdready.pdb --keep-ligands
    python scripts/md_ready.py rec.pdb --dry-run
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_protein import ALIASES, HEAVY_ATOMS, WATER  # noqa: E402

# 修饰残基 → 母体标准残基
TO_STANDARD = {
    "MSE": "MET", "FME": "MET", "TPO": "THR", "SEP": "SER", "PTR": "TYR",
    "CSO": "CYS", "CSD": "CYS", "CME": "CYS", "OCY": "CYS", "CSS": "CYS",
    "CYX": "CYS", "CYM": "CYS", "KCX": "LYS", "LLP": "LYS", "MLY": "LYS",
    "M3L": "LYS", "ALY": "LYS", "LYN": "LYS", "HYP": "PRO", "PCA": "GLU",
    "OAS": "SER", "OAT": "ALA", "NLE": "LEU", "SAR": "GLY", "ABA": "ALA",
    "DAL": "ALA", "DLE": "LEU", "DVA": "VAL", "ASH": "ASP", "GLH": "GLU",
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS",
    "HSP": "HIS",
}
STANDARD = set(HEAVY_ATOMS)
BACKBONE = {"N", "CA", "C", "O", "OXT", "OT1", "OT2"}


def _say(*a):
    print(*a, flush=True)


def residues(path: Path):
    """按出现顺序给出 [(chain, resseq+icode, resname, [原子行])]。

    同一残基里同名的原子只留第一个：修饰残基常同时以 ATOM 和 HETATM 各写一遍
    （AKT1 的 TPO308、PTGS2 的 OAS530 都是），重复原子坐标完全重合，进 mdrun 直接
    把最大力推到 inf，最小化当场崩掉。
    """
    out: list[list] = []
    for ln in path.read_text(errors="replace").splitlines():
        if not ln.startswith(("ATOM", "HETATM")) or len(ln) < 27:
            continue
        key = (ln[21], ln[22:27])
        if not out or out[-1][0] != key:
            out.append([key, ln[17:20].strip(), [], set()])
        seen = out[-1][3]
        name = ln[12:16].strip()
        if name in seen:
            continue
        seen.add(name)
        out[-1][2].append(ln)
    return [(k, n, a) for k, n, a, _ in out]


def _renamed_atoms(atom_lines: list[str], std: str) -> tuple[list[str], int]:
    allowed = set(HEAVY_ATOMS.get(std, [])) | BACKBONE
    alias = ALIASES.get(std, (std, {}))[1]      # 例如 MSE 的 SE→SD
    keep, dropped = [], 0
    for ln in atom_lines:
        name = ln[12:16].strip()
        if name in alias:
            name = alias[name]
            ln = f"{ln[:12]}{name:>4}{ln[16:]}"
        if name in allowed:
            keep.append(ln)
        else:
            dropped += 1
    return keep, dropped


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    ap = argparse.ArgumentParser(description="受体结构 MD 前清理：修饰残基降级、去小分子")
    ap.add_argument("pdb")
    ap.add_argument("--out", help="输出 PDB；默认 <输入>_mdready.pdb")
    ap.add_argument("--keep-ligands", action="store_true", help="保留共晶小分子与离子")
    ap.add_argument("--keep-water", action="store_true", help="保留水")
    ap.add_argument("--dry-run", action="store_true", help="只报告要改什么")
    a = ap.parse_args(argv)

    src = Path(a.pdb)
    if not src.exists():
        sys.exit(f"找不到 {src}")
    out = Path(a.out) if a.out else src.with_name(src.stem + "_mdready.pdb")

    res = residues(src)
    # 每条链最后一个残基的下标（OXT 只允许出现在这里）
    last_of_chain: dict[str, int] = {}
    for i, ((chain, _key), name, _atoms) in enumerate(res):
        if name not in WATER:
            last_of_chain[chain] = i

    kept_lines, dropped_res, renamed = [], [], []
    for ri, ((chain, key), name, atoms) in enumerate(res):
        if ri != last_of_chain.get(chain):
            # MODELLER 会把模板 C 端的 OXT 抄到每段末尾；不是链尾就不能留，
            # 否则 pdb2gmx 报 "Atom OXT ... was not found in rtp entry"
            atoms = [ln for ln in atoms if ln[12:16].strip() not in ("OXT", "OT1", "OT2")]
        # 要有 N、CA、C 三件才算肽链残基。小分子常带叫 N、C 的原子（STAT3 的
        # KQV 就是），只查"有没有 N 或 C"会把它当成蛋白残基留下来
        names = {ln[12:16].strip() for ln in atoms}
        has_backbone = {"N", "CA", "C"} <= names
        if name in WATER:
            if a.keep_water:
                kept_lines += atoms
            else:
                dropped_res.append(f"{name} {chain}{key.strip()}(水)")
            continue
        if name in STANDARD:
            kept_lines += atoms
            continue
        if name in TO_STANDARD:
            std = TO_STANDARD[name]
            keep, n_drop = _renamed_atoms(atoms, std)
            renamed.append(f"{name}{key.strip()} → {std}"
                           + (f"（删 {n_drop} 个原子）" if n_drop else ""))
            kept_lines += [f"{ln[:17]}{std:>3}{ln[20:]}" for ln in keep]
            continue
        # 其余：小分子、离子、未知残基
        if a.keep_ligands:
            kept_lines += atoms
        elif has_backbone:
            # 有主链却不在标准表里——多半是没登记的新修饰残基，别悄悄删
            kept_lines += atoms
            dropped_res.append(f"! {name} {chain}{key.strip()} 有主链但不认识，原样保留")
        else:
            dropped_res.append(f"{name} {chain}{key.strip()}")

    _say(f"{src.name}：{len(res)} 个残基")
    if renamed:
        _say(f"  修饰残基降级 {len(renamed)} 处：")
        for r in renamed[:12]:
            _say(f"    {r}")
        if len(renamed) > 12:
            _say(f"    ……共 {len(renamed)} 处")
    if dropped_res:
        _say(f"  删除 {len(dropped_res)} 个残基：")
        for r in dropped_res[:12]:
            _say(f"    {r}")
        if len(dropped_res) > 12:
            _say(f"    ……共 {len(dropped_res)} 个")
    if not renamed and not dropped_res:
        _say("  没什么要改的，原样输出")

    if a.dry_run:
        _say("（--dry-run，没有写文件）")
        return 0

    head = [ln for ln in src.read_text(errors="replace").splitlines()
            if not ln.startswith(("ATOM", "HETATM", "TER", "END", "MASTER"))]
    # 链之间补 TER：删掉小分子后原来的 TER 位置没了，靠 chain ID 变化重插
    body, prev_chain = [], None
    for ln in kept_lines:
        ch = ln[21] if len(ln) > 21 else " "
        if prev_chain is not None and ch != prev_chain:
            body.append("TER")
        body.append(ln)
        prev_chain = ch
    body.append("TER")
    out.write_text("\n".join(head + body + ["END"]) + "\n", encoding="utf-8")
    _say(f"写出 {out}")
    _say(f'复查：python scripts/check_protein.py "{out}"')
    return 0


if __name__ == "__main__":
    sys.exit(main())
