#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
split_complex.py —— 把"蛋白 + 配体"的复合物 PDB 拆成两份输入文件

有些对接网站（CB-Dock、SwissDock 等）只给一个复合物 PDB。MD 流水线要的是
分开的两份：蛋白给 pdb2gmx，配体给 GAFF2 参数化。拆分不移动任何原子。

命令行用法
----------
    python scripts/split_complex.py --in complex.pdb --out-dir split

    # 复合物里有多个小分子时，指定要哪一个
    python scripts/split_complex.py --in complex.pdb --out-dir split --lig-resname UNL

输出：split/protein.pdb（蛋白原子 + SSBOND 表头）
      split/ligand.mol2（补全氢，可直接作为流水线的"配体位姿"）

也作为库被 mdpipe.py 调用（--complex 与 inbox 里的复合物自动拆分）。
"""
import argparse
import subprocess
import sys
from collections import Counter
from pathlib import Path

WATER = {"HOH", "WAT", "DOD", "SOL"}

# 标准氨基酸 + 常见变体；不在这里的 ATOM 残基八成是配体
STD_RES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN", "ACE", "NME", "UNK",
    "MSE", "SEP", "TPO", "PTR", "CSO", "KCX", "MLY", "PCA", "FME", "HYP",
}


def find_ligand(lines, resname: str | None = None) -> tuple[str, bool]:
    """找出配体的残基名。返回 (残基名, 配体是否以 ATOM 记录写出)。

    优先非水 HETATM；没有 HETATM 时，看有没有唯一的非标准 ATOM 残基
    （有些软件把配体写成 ATOM，残基名是 *** 之类）。找不到抛 ValueError。
    """
    het = Counter()
    for l in lines:
        if l.startswith("HETATM"):
            name = l[17:20].strip().upper()
            if name not in WATER:
                het[name] += 1
    odd = Counter(l[17:20].strip().upper() for l in lines
                  if l.startswith("ATOM") and l[17:20].strip().upper() not in STD_RES)

    if resname:
        lig = resname.strip().upper()
        allres = Counter(l[17:20].strip().upper() for l in lines
                         if l.startswith(("ATOM", "HETATM")))
        if lig not in allres:
            raise ValueError(f"没有找到残基 {lig}；非水 HETATM：{dict(het)}，"
                             f"非标准 ATOM 残基：{dict(odd)}")
    elif het and len(het) == 1:
        lig = next(iter(het))
    elif odd and len(odd) == 1:
        lig = next(iter(odd))
    elif het or odd:
        raise ValueError(f"复合物里有多个候选小分子，用 --lig-resname 指定一个："
                         f"{ {**dict(het), **dict(odd)} }")
    else:
        raise ValueError("没有找到配体：既没有非水 HETATM，也没有非标准 ATOM 残基。")
    return lig, lig in odd


def split_lines(lines, lig: str, want_atom: bool):
    """切成 (蛋白行, 配体行)。want_atom 表示配体写在 ATOM 记录里。"""
    prot = [l for l in lines
            if l.startswith(("ATOM", "SSBOND"))
            and not (want_atom and l.startswith("ATOM")
                     and l[17:20].strip().upper() == lig)]
    lig_lines = [l for l in lines
                 if (l.startswith("HETATM") or (want_atom and l.startswith("ATOM")))
                 and l[17:20].strip().upper() == lig]
    return prot, lig_lines


def split_to(src, out_dir, resname: str | None = None) -> tuple[Path, Path]:
    """拆分 + 给配体补氢。返回 (protein.pdb, ligand.mol2)。"""
    src, out = Path(src), Path(out_dir)
    lines = src.read_text(encoding="utf-8", errors="replace").splitlines()
    lig, want_atom = find_ligand(lines, resname)
    prot, lig_lines = split_lines(lines, lig, want_atom)
    if not prot:
        raise RuntimeError(f"{src} 里没有蛋白原子")

    out.mkdir(parents=True, exist_ok=True)
    prot_path = out / "protein.pdb"
    prot_path.write_text("\n".join(prot + ["TER", "END"]) + "\n", encoding="utf-8")

    tmp = out / "_ligand_raw.pdb"
    body = [f"{l[:6]}{i:5d}{l[11:17]}LIG{l[20:]}" for i, l in enumerate(lig_lines, 1)]
    tmp.write_text("\n".join(body + ["END"]) + "\n", encoding="utf-8")

    mol2 = out / "ligand.mol2"
    try:
        r = subprocess.run(["obabel", str(tmp), "-O", str(mol2), "-h"],
                           capture_output=True, text=True)
    finally:
        tmp.unlink(missing_ok=True)
    if r.returncode != 0 or not mol2.exists():
        raise RuntimeError(f"obabel 补氢失败：\n{r.stderr}")
    return prot_path, mol2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="复合物 PDB")
    ap.add_argument("--out-dir", default="split", help="输出目录")
    ap.add_argument("--lig-resname", default=None,
                    help="配体的残基名；不给则在只有一种小分子时自动选")
    a = ap.parse_args()

    src = Path(a.inp)
    lines = src.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        lig, want_atom = find_ligand(lines, a.lig_resname)
    except ValueError as e:
        sys.exit(str(e))
    if want_atom:
        print(f"配体以 ATOM 记录给出，残基名 {lig}，按它处理")

    try:
        prot_path, mol2 = split_to(src, a.out_dir, a.lig_resname)
    except RuntimeError as e:
        sys.exit(str(e))

    n_prot = sum(1 for l in prot_path.read_text(encoding="utf-8").splitlines()
                 if l.startswith("ATOM"))
    n_lig = sum(1 for l in lines
                if l[17:20].strip().upper() == lig
                and (l.startswith("HETATM") or (want_atom and l.startswith("ATOM"))))
    print(f"蛋白 -> {prot_path}（{n_prot} 个原子）")
    print(f"配体 -> {mol2}（残基 {lig}，{n_lig} 个重原子，已补氢）")
    print("下一步：把这两个文件分别作为「蛋白 PDB」和「配体位姿」投进流水线。")


if __name__ == "__main__":
    main()
