#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""蛋白结构体检：缺残基、断链、altLoc、非标准残基、缺重原子、杂原子。

为什么需要这一步：对接下载的受体 PDB 常常缺一段无序区（REMARK 465）。
直接丢进 pdb2gmx，缺口处会被当成新链的端点补上带电末端，蛋白内部凭空多出
一对 NH3+/COO-，跑起来这一段会散开。跑之前先体检，缺口长的补链、短的标出来。

用法：
    python scripts/check_protein.py rec.pdb                 # 打印摘要
    python scripts/check_protein.py rec.pdb --md rec_体检.md  # 同时写报告
    python scripts/check_protein.py *.pdb --json out.json   # 批量，出 JSON

零依赖（只用标准库）。.cif 会先用 OpenBabel 转成 PDB 再读。
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# ---------------------------------------------------------------- 残基定义

# 标准氨基酸的重原子清单（不含氢）
HEAVY_ATOMS = {
    "ALA": "N CA C O CB",
    "ARG": "N CA C O CB CG CD NE CZ NH1 NH2",
    "ASN": "N CA C O CB CG OD1 ND2",
    "ASP": "N CA C O CB CG OD1 OD2",
    "CYS": "N CA C O CB SG",
    "GLN": "N CA C O CB CG CD OE1 NE2",
    "GLU": "N CA C O CB CG CD OE1 OE2",
    "GLY": "N CA C O",
    "HIS": "N CA C O CB CG ND1 CD2 CE1 NE2",
    "ILE": "N CA C O CB CG1 CG2 CD1",
    "LEU": "N CA C O CB CG CD1 CD2",
    "LYS": "N CA C O CB CG CD CE NZ",
    "MET": "N CA C O CB CG SD CE",
    "PHE": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ",
    "PRO": "N CA C O CB CG CD",
    "SER": "N CA C O CB OG",
    "THR": "N CA C O CB OG1 CG2",
    "TRP": "N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH",
    "VAL": "N CA C O CB CG1 CG2",
}
HEAVY_ATOMS = {k: set(v.split()) for k, v in HEAVY_ATOMS.items()}

# 常见变体名 -> (标准名, 需要替换的原子名表)
ALIASES = {
    "HID": ("HIS", {}), "HIE": ("HIS", {}), "HIP": ("HIS", {}), "HSD": ("HIS", {}),
    "HSE": ("HIS", {}), "HSP": ("HIS", {}),
    "CYX": ("CYS", {}), "CYM": ("CYS", {}),
    "ASH": ("ASP", {}), "GLH": ("GLU", {}), "LYN": ("LYS", {}),
    "MSE": ("MET", {"SE": "SD"}),
    "SEC": ("CYS", {"SE": "SG"}),
    "UNK": ("UNK", {}),
}
STANDARD = set(HEAVY_ATOMS)

WATER = {"HOH", "WAT", "DOD", "H2O", "TIP", "SOL"}
COMMON_IONS = {"NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE", "CU", "NI", "CD",
               "SO4", "PO4", "IOD", "BR"}

# Amber/GROMACS 的原子名和 IUPAC 不完全一样，比原子齐全度时要归一。
# 最典型的是异亮氨酸：IUPAC 叫 CD1，Amber（ff14SB/ff19SB）叫 CD。
AMBER_ALIAS = {"ILE": {"CD": "CD1"}}

BACKBONE = {"N", "CA", "C", "O"}
PEPTIDE_BOND_MAX = 1.7   # C(i)-N(i+1) 超过这个距离算断链（正常 1.33 Å）


def canon(name: str) -> str:
    name = name.strip().upper()
    return ALIASES.get(name, (name, {}))[0]


# ---------------------------------------------------------------- 读结构

def cif_to_pdb(path: Path) -> Path:
    """用 OpenBabel 把 mmCIF 转成 PDB（临时文件）。"""
    obabel = shutil.which("obabel")
    if not obabel:
        sys.exit("读 .cif 需要 OpenBabel（obabel），PATH 里没找到；先用其它工具转成 .pdb")
    out = Path(tempfile.mkdtemp()) / (path.stem + ".pdb")
    r = subprocess.run([obabel, str(path), "-O", str(out)], capture_output=True, text=True)
    if not out.exists():
        sys.exit(f"OpenBabel 转换失败：{path}\n{r.stderr[-500:]}")
    return out


def parse_pdb(path: Path) -> dict:
    """解析 PDB 的 ATOM/HETATM/SEQRES/REMARK 465。"""
    chains: dict[str, list[dict]] = {}      # 链 -> 残基列表（文件顺序）
    seqres: dict[str, list[tuple[str, int]]] = {}
    remark465: list[tuple[str, int, str]] = []
    het: dict[str, int] = {}
    het_res: list[tuple[str, int, str]] = []     # (链, 残基号, 残基名)，用来识别链内修饰残基
    n_altloc = 0
    n_lowocc = 0
    has_h = False
    model_count = 0

    with path.open("r", errors="replace") as fh:
        for line in fh:
            rec = line[:6].strip()

            if rec == "MODEL":
                model_count += 1
                if model_count > 1:          # 只取第一个模型
                    break
                continue

            if rec == "SEQRES":
                ch = line[11]
                for i in range(13, 13 + 4 * 13, 4):
                    name = line[i:i + 3].strip() if i + 3 <= len(line) else ""
                    if not name:
                        break
                    seqres.setdefault(ch, []).append(name)
                continue

            if line.startswith("REMARK 465"):
                # REMARK 465   M RES C SSSEQI
                ch = line[19] if len(line) > 19 else " "
                try:
                    num = int(line[21:26])
                    name = line[15:18].strip()
                except (ValueError, IndexError):
                    continue
                remark465.append((ch, num, name))
                continue

            if rec not in ("ATOM", "HETATM"):
                continue

            name = line[12:16].strip()
            alt = line[16]
            resn = line[17:20].strip().upper()
            ch = line[21]
            icode = line[26] if len(line) > 26 else " "

            if name.startswith("H") or (len(name) > 1 and name[0].isdigit() and name[1] == "H"):
                has_h = True

            if resn in WATER:
                het["水"] = het.get("水", 0) + 1
                continue

            try:
                rs0 = int(line[22:26])
            except ValueError:
                rs0 = None

            is_polymer = rec == "ATOM" or canon(resn) in STANDARD
            if not is_polymer:
                het[resn] = het.get(resn, 0) + 1
                if rec == "HETATM" and rs0 is not None and (ch, rs0, resn) not in het_res:
                    het_res.append((ch, rs0, resn))
                continue

            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                occ = float(line[54:60]) if line[54:60].strip() else 1.0
            except ValueError:
                continue

            if alt and alt != " ":
                n_altloc += 1
            if occ < 0.5:
                n_lowocc += 1

            rs = int(line[22:26])
            lst = chains.setdefault(ch, [])
            # 同一 (残基号, 插入码) 归为一个残基
            if lst and lst[-1]["num"] == rs and lst[-1]["icode"] == icode:
                r = lst[-1]
            else:
                r = {"num": rs, "icode": icode, "name": resn, "atoms": {}}
                lst.append(r)
            # altLoc 只保留 A（或第一个）
            if alt in (" ", "A") and name not in r["atoms"]:
                r["atoms"][name] = (x, y, z)

    return {"chains": chains, "seqres": seqres, "remark465": remark465,
            "het": het, "het_res": het_res,
            "n_altloc": n_altloc, "n_lowocc": n_lowocc, "has_h": has_h}


# ---------------------------------------------------------------- 检查

def analyse(path: Path) -> dict:
    src = path
    if path.suffix.lower() in (".cif", ".mmcif"):
        src = cif_to_pdb(path)

    d = parse_pdb(src)
    chains = d["chains"]
    issues: list[dict] = []
    chain_reports = []

    for ch, res in chains.items():
        if not res:
            continue
        nums_present = {r["num"] for r in res}

        # 链内修饰残基：HETATM 记的（TPO/SEP/MLY…），编号正好插在现有残基中间。
        # 它们占一个序列位置，所以后面拼"序列"时要算进去，否则会误报成缺残基。
        chain_het = [(n, nm) for (c, n, nm) in d["het_res"]
                     if c == ch and n not in nums_present
                     and ((n - 1) in nums_present or (n + 1) in nums_present)]

        seq_items = ([("std", r["num"], r) for r in res]
                     + [("het", n, nm) for n, nm in chain_het])
        seq_items.sort(key=lambda t: t[2]["num"] if t[0] == "std" else t[1])

        # 1) 编号断口（het 残基只占位，不参与断口判断）
        gaps_num = []
        for (_, na, _), (_, nb, _) in zip(seq_items, seq_items[1:]):
            if nb - na > 1:
                gaps_num.append({"before": na, "after": nb, "n": nb - na - 1})

        # 2) 几何断口：只有相邻两位都有坐标时才比 C–N
        breaks_geo = []
        for (ka, na, ra), (kb, nb, rb) in zip(seq_items, seq_items[1:]):
            if ka != "std" or kb != "std":
                continue
            ca, n_atom = ra["atoms"].get("C"), rb["atoms"].get("N")
            if ca and n_atom:
                dist = math.dist(ca, n_atom)
                if dist > PEPTIDE_BOND_MAX:
                    breaks_geo.append({"before": na, "after": nb, "dist": round(dist, 2),
                                       "consecutive": nb - na == 1})

        # 3) 非标准残基与缺重原子
        incomplete = []
        nonstd = [f'{r["name"]}{r["num"]}' for r in res
                  if canon(r["name"]) not in STANDARD]
        for r in res:
            base = canon(r["name"])
            if base not in HEAVY_ATOMS:
                continue
            mapping = dict(ALIASES.get(r["name"].upper(), ("", {}))[1])
            mapping.update(AMBER_ALIAS.get(base, {}))
            expected = {mapping.get(a, a) for a in HEAVY_ATOMS[base]}
            have = {mapping.get(a, a) for a in r["atoms"]}
            missing = expected - have
            # 末端残基允许缺一个氧（OXT 不是必需），侧链缺失才记
            missing -= {"OXT"}
            if missing and missing != {"O"}:
                incomplete.append({"resid": f'{r["name"]}{r["num"]}',
                                   "missing": sorted(missing)})

        remark = [(n, nm) for (c, n, nm) in d["remark465"] if c == ch]

        chain_reports.append({
            "chain": ch,
            "residues": len(seq_items),
            "range": f'{res[0]["num"]}-{res[-1]["num"]}',
            "seq": "".join(three_to_one(canon(r["name"])) for r in res),
            "gaps_numbering": gaps_num,
            "breaks_geometry": breaks_geo,
            "missing_in_remark465": len(remark),
            "nonstandard": (nonstd + [f"{nm}{n}" for n, nm in chain_het])[:20],
            "incomplete_residues": incomplete[:40],
        })

        for g in gaps_num:
            issues.append({"level": "ERROR" if g["n"] > 1 else "WARN", "chain": ch,
                           "what": f'缺 {g["n"]} 个残基（{g["before"]}–{g["after"]} 之间）'})
        for b in breaks_geo:
            if b["consecutive"]:
                # 编号连续但键长偏长：多半是补链/建模工具留下的接缝。
                # pdb2gmx 按残基序列建键、不看坐标，拓扑是连通的，EM 会把它收拢，
                # 所以不算致命问题。
                issues.append({"level": "WARN", "chain": ch,
                               "what": f'接缝键长偏长 {b["before"]}–{b["after"]}'
                                       f'（C–N {b["dist"]} Å，EM 会收拢）'})
            else:
                issues.append({"level": "ERROR", "chain": ch,
                               "what": f'链在 {b["before"]}–{b["after"]} 之间断开'
                                       f'（C–N {b["dist"]} Å）'})
        for it in incomplete:
            issues.append({"level": "WARN", "chain": ch,
                           "what": f'{it["resid"]} 缺重原子 {",".join(it["missing"])}'})
        for n in nonstd:
            issues.append({"level": "WARN", "chain": ch, "what": f"非标准残基 {n}"})

        # 链内修饰残基（TPO/SEP/MLY…）在链中间占着一整个残基位置，
        # pdb2gmx 认不出来会直接报 Residue 'TPO' not found
        for hnum, hname in chain_het:
            issues.append({"level": "ERROR", "chain": ch,
                           "what": f"链内修饰残基 {hname}{hnum}（pdb2gmx 不认，"
                                   f"需先换成标准残基或单独参数化）"})

    if d["n_altloc"]:
        issues.append({"level": "WARN", "chain": "*",
                       "what": f'{d["n_altloc"]} 个原子有 altLoc（已只取 A）'})
    if d["n_lowocc"]:
        issues.append({"level": "INFO", "chain": "*",
                       "what": f'{d["n_lowocc"]} 个原子占有率 < 0.5'})

    # 判定
    n_err = sum(1 for i in issues if i["level"] == "ERROR")
    n_warn = sum(1 for i in issues if i["level"] == "WARN")
    n_seam = sum(1 for i in issues if str(i["what"]).startswith("接缝键长"))

    if n_err:
        verdict, advice = "FAIL", "有断链或整段缺残基，直接 pdb2gmx 会在缺口处插带电末端。先补链（complete_protein.py），或只保留完整的那条链。"
    elif n_warn and n_seam == n_warn:
        verdict, advice = "WARN", "只剩补链接缝的键长偏长。pdb2gmx 按残基序列建键、不看坐标，拓扑是连通的，跑 EM 就收拢了，可以直接进流水线。"
    elif n_warn:
        verdict, advice = "WARN", "可以跑，但建议先补齐侧链重原子；补链工具顺带能做。"
    else:
        verdict, advice = "PASS", "结构完整，可直接进流水线。"

    return {"file": str(path), "verdict": verdict, "advice": advice,
            "n_chains": len([c for c in chain_reports]),
            "chains": chain_reports, "issues": issues,
            "hetero": d["het"], "has_hydrogens": d["has_h"],
            "remark465_total": len(d["remark465"]), "seqres_chains": list(d["seqres"])}


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "UNK": "X",
}


def three_to_one(n: str) -> str:
    return THREE_TO_ONE.get(n, "X")


# ---------------------------------------------------------------- 输出

def to_markdown(reports: list[dict]) -> str:
    out = ["# 蛋白结构体检报告", ""]
    for r in reports:
        out.append(f'## {Path(r["file"]).name} — **{r["verdict"]}**')
        out.append("")
        out.append(r["advice"])
        out.append("")
        out.append("| 链 | 残基数 | 编号范围 | 编号缺口 | 几何断链 | REMARK465 缺残基 |")
        out.append("|---|---|---|---|---|---|")
        for c in r["chains"]:
            gaps = "、".join(f'{g["before"]}→{g["after"]}(缺{g["n"]})'
                            for g in c["gaps_numbering"]) or "—"
            brk = "、".join(f'{b["before"]}→{b["after"]}' for b in c["breaks_geometry"]) or "—"
            out.append(f'| {c["chain"]} | {c["residues"]} | {c["range"]} | {gaps} | {brk} | '
                       f'{c["missing_in_remark465"]} |')
        out.append("")
        if r["issues"]:
            out.append("| 级别 | 链 | 内容 |")
            out.append("|---|---|---|")
            for i in r["issues"][:60]:
                out.append(f'| {i["level"]} | {i["chain"]} | {i["what"]} |')
            if len(r["issues"]) > 60:
                out.append(f'| … | | 另有 {len(r["issues"]) - 60} 条 |')
            out.append("")
        if r["hetero"]:
            het = "、".join(f"{k}×{v}" for k, v in sorted(r["hetero"].items(), key=lambda x: -x[1])[:12])
            out.append(f"非蛋白组分：{het}")
            out.append("")
    return "\n".join(out)


def to_text(r: dict) -> str:
    lines = [f'【{Path(r["file"]).name}】{r["verdict"]}  {r["advice"]}']
    for c in r["chains"]:
        gap = ""
        if c["gaps_numbering"]:
            gap = " 缺 " + "、".join(f'{g["before"]}→{g["after"]}({g["n"]})'
                                     for g in c["gaps_numbering"])
        brk = ""
        if c["breaks_geometry"]:
            brk = " 断链@" + "、".join(f'{b["before"]}→{b["after"]}' for b in c["breaks_geometry"])
        inc = f' 缺重原子×{len(c["incomplete_residues"])}' if c["incomplete_residues"] else ""
        lines.append(f'  链 {c["chain"]}: {c["residues"]} 残基 ({c["range"]}){gap}{brk}{inc}')
    errs = [i for i in r["issues"] if i["level"] == "ERROR"]
    for i in errs[:8]:
        lines.append(f'  ! {i["chain"]}: {i["what"]}')
    if len(errs) > 8:
        lines.append(f'  ! 另有 {len(errs) - 8} 条错误')
    return "\n".join(lines)


def main(argv=None):
    # Windows 控制台默认 GBK，遇到 Å 之类字符会崩；统一走 UTF-8
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    ap = argparse.ArgumentParser(
        description="蛋白结构体检：缺残基、断链、altLoc、非标准残基、缺重原子")
    ap.add_argument("files", nargs="+", help="PDB / mmCIF 文件")
    ap.add_argument("--md", help="同时写 Markdown 报告到这个路径")
    ap.add_argument("--json", help="同时写 JSON 到这个路径")
    a = ap.parse_args(argv)

    reports = []
    for f in a.files:
        p = Path(f)
        if not p.exists():
            print(f"跳过（不存在）：{f}", file=sys.stderr)
            continue
        try:
            reports.append(analyse(p))
        except Exception as e:                     # 单个文件坏了不影响批量
            print(f"跳过（解析失败）：{f}  {e}", file=sys.stderr)

    for r in reports:
        print(to_text(r))
        print()

    n_fail = sum(1 for r in reports if r["verdict"] == "FAIL")
    n_warn = sum(1 for r in reports if r["verdict"] == "WARN")
    print(f"共 {len(reports)} 个文件：{n_fail} 个 FAIL，{n_warn} 个 WARN，"
          f"{len(reports) - n_fail - n_warn} 个 PASS")

    if a.md:
        Path(a.md).write_text(to_markdown(reports), encoding="utf-8")
        print(f"报告已写入 {a.md}")
    if a.json:
        Path(a.json).write_text(json.dumps(reports, ensure_ascii=False, indent=2),
                                encoding="utf-8")
        print(f"JSON 已写入 {a.json}")

    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
