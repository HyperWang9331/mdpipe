#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 MODELLER 补缺失的 loop：晶体结构原位不动，只对缺口建模。

为什么还要有这个——PDBFixer 只能从 SEQRES 推缺口，遇到编号重复、SEQRES 被
清洗掉、或者 SEQRES 本身不全的文件就完全失灵。MODELLER 是"模板建模"，给出的
loop 构象质量也更好，代价是慢（每秒几十个残基）和需要 license。

做法（把网上手工教程整套自动化）：
  1. 把目标链单独抠出来，用它的已有残基当模板（structureX）。
  2. 用 SEQRES 给出完整序列（sequence），两行做 PIR 比对，缺口位置填 '-'。
  3. select_atoms 只选缺口那几段 → 已有残基一个原子都不动。
  4. 跑 N 个模型，按 molpdf 选最好的。
  5. 多链结构把补好的链拼回原来的结构里，残基编号还原。

晶体残基 ↔ SEQRES 的对位有两套办法，先试编号、不行再按片段找：
  * 编号法：SEQRES 下标与残基号的差投票，适合 STAT3/AKT1 这种编号规整的文件。
  * 片段法：用主链 C–N 距离和编号跳变找出真实断口，把链切成几段，每段单独在
    SEQRES 里滑窗找位。BCL2 的 39→89 这种编号和序列完全脱节的文件只能这么办。
纯序列 LCS 是不能用的：缺一个 ASN 时它会挑一个同名的 ASN 去配，把缺口错切到
相邻残基中间，补出来的 loop 几何上根本不成立（做成 43 Å 的断链）。

MODELLER 内部一律按"位置号"（目标序列第几个残基，1 起）工作，模板 PDB 也照位置
号重编，避免原文件编号乱跳导致 structureX 对不上；输出时再换算回原编号。

license：MODELLER 学术免费，但要注册领 key（https://salilab.org/modeller/）。
本流水线不附带 MODELLER，也不代填 license；使用前按官方说明把 key 写进：
    <env>/Library/modeller/modlib/modeller/config.py 里 license = r'你的KEY'

用法：
    python scripts/fill_loops_modeller.py rec.pdb --out rec_modeled.pdb
    python scripts/fill_loops_modeller.py rec.pdb --out rec_modeled.pdb --models 10
    python scripts/fill_loops_modeller.py rec.pdb --chain A --models 3
    python scripts/fill_loops_modeller.py rec.pdb --build-termini     # 首尾悬空段也建

资源：单线程 + BelowNormal，不碰 GPU，可以和正在跑的 MD 并行。
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from complete_protein import _lcs_pairs, _seqres_by_chain  # noqa: E402
from check_protein import HEAVY_ATOMS  # noqa: E402
from md_ready import TO_STANDARD  # noqa: E402

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

MODIFIED_TO_STD = {
    "MSE": "MET", "FME": "MET", "SEP": "SER", "TPO": "THR", "PTR": "TYR",
    "CSO": "CYS", "CSD": "CYS", "CME": "CYS", "OCY": "CYS", "CSS": "CYS",
    "KCX": "LYS", "LLP": "LYS", "MLY": "LYS", "M3L": "LYS", "ALY": "LYS",
    "HYP": "PRO", "PCA": "GLU", "OAS": "SER", "OAT": "ALA", "NLE": "LEU",
    "SAR": "GLY", "ABA": "ALA", "DAL": "ALA", "DLE": "LEU", "DVA": "VAL",
}

PEPTIDE_MAX = 2.0   # C–N 超过这个距离就认定是真实断口


def _say(*a):
    print(*a, flush=True)


def find_modeller_python() -> Path | None:
    for base in (Path.home() / "miniconda3" / "envs", Path.home() / "anaconda3" / "envs"):
        p = base / "modeller" / "python.exe"
        if p.exists():
            try:
                r = subprocess.run([str(p), "-c", "import modeller"],
                                   capture_output=True, timeout=180)
                if r.returncode == 0:
                    return p
            except Exception:
                pass
    return None


def read_chain(pdb: Path, chain: str) -> tuple[list[str], list[tuple[int, str]], bool]:
    """取出某条链的 ATOM 行，以及按顺序的 (残基号, 残基名)。

    有插入码（70A 这种）的链，顺着顺序重编成连续号：MODELLER 只认整数残基号，
    70 和 70A 在它看来是同一个残基，不改就会报 Sequence difference。
    """
    raw = [ln for ln in pdb.read_text(errors="replace").splitlines()
           if ln.startswith(("ATOM  ", "HETATM")) and ln[21] == chain]
    groups: dict[tuple[int, str], list[str]] = {}
    order: list[tuple[int, str]] = []
    for ln in raw:
        try:
            num = int(ln[22:26])
        except ValueError:
            continue
        key = (num, ln[26] if len(ln) > 26 else " ")
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(ln)

    # 只留链上的残基（有 N/CA/C 三件）。修饰残基常以 HETATM 写（TPO、OAS），
    # 只读 ATOM 会把它当缺口让 MODELLER 重建一遍，输出里就和原残基叠在一起，
    # 能量最小化直接 inf。这里统一降级成母体标准残基后再参与比对。
    kept, unknown = [], []
    for key in order:
        grp, seen_atoms = [], set()
        for x in groups[key]:
            an = x[12:16].strip()
            if an in seen_atoms:
                continue
            seen_atoms.add(an)
            grp.append(x)
        names = {x[12:16].strip() for x in grp}
        if not {"N", "CA", "C"} <= names:
            continue                      # 水、离子、小分子，不进链
        resname = grp[0][17:20].strip().upper()
        lines = grp
        if resname in HEAVY_ATOMS:
            std = resname
        else:
            std = TO_STANDARD.get(resname)
            if not std:
                unknown.append(resname)
                std = resname
            allowed = set(HEAVY_ATOMS.get(std, ())) | {"N", "CA", "C", "O",
                                                       "OXT", "OT1", "OT2"}
            lines = [f"ATOM  {x[6:17]}{std:>3}{x[20:]}" for x in grp
                     if x[12:16].strip() in allowed]
        kept.append((key, std, lines))

    if not kept:
        return [], [], False
    if unknown:
        _say(f"  链 {chain}：不认识的残基 {sorted(set(unknown))}，原样参与建模")

    has_icode = any(k[1].strip() for k, _, _ in kept)
    base = kept[0][0][0]
    if not has_icode:
        return ([x for _, _, ls in kept for x in ls],
                [(k[0], s) for k, s, _ in kept], False)

    # 有插入码：按出现顺序重编成连续号（MODELLER 只认整数残基号）
    out = []
    for i, (_k, _s, ls) in enumerate(kept):
        out += [f"{x[:22]}{base + i:4d} {x[27:]}" for x in ls]
    return out, [(base + i, s) for i, (_k, s, _l) in enumerate(kept)], True


def _same_res(seqres_name: str, obs_name: str) -> bool:
    return seqres_name == MODIFIED_TO_STD.get(obs_name.upper(), obs_name.upper())


def _aa1(name: str) -> str | None:
    """三字母 → 单字母。PIR 的序列行必须是单字母，写三字母 MODELLER 会把首字母
    当成一个"未知残基"来读（报 Unknown residue type ... U）。"""
    return THREE_TO_ONE.get(MODIFIED_TO_STD.get(name.upper(), name.upper()))


def _iter_fasta(path: Path):
    head, buf = None, []
    for ln in path.read_text(errors="replace").splitlines():
        if ln.startswith(">"):
            if head is not None:
                yield head, "".join(buf)
            head, buf = ln[1:].strip(), []
        elif head is not None:
            buf.append(ln.strip())
    if head is not None:
        yield head, "".join(buf)


def _load_fasta(paths) -> tuple[dict[str, list[str]], list[str] | None]:
    """读 FASTA，给出 (按链分配的序列, 默认序列)。

    记录名是单个字母（A/B/C）就分配到那条链；否则当作默认序列，给所有没有
    SEQRES 的链共用——PyMOL/AutoDock 清洗过的受体常常只剩一条链，也没了
    SEQRES，这时一条 UniProt 序列就够。
    """
    to3 = {v: k for k, v in THREE_TO_ONE.items()}
    by_chain: dict[str, list[str]] = {}
    default: list[str] | None = None
    for p in paths:
        p = Path(p)
        if not p.exists():
            sys.exit(f"找不到序列文件 {p}")
        for head, seq in _iter_fasta(p):
            residues = []
            for ch in seq.upper():
                if ch in to3:
                    residues.append(to3[ch])
                elif ch in "-.":
                    continue
                else:
                    sys.exit(f"{p.name}：序列里有不认识的字符 {ch!r}（只能有 20 种标准氨基酸）")
            tag = head.split()[0] if head.split() else ""
            if len(tag) == 1 and tag.isalpha():
                by_chain[tag.upper()] = residues
            else:
                default = residues
    return by_chain, default


def _break_positions(lines: list[str]) -> list[int]:
    """找出 file 顺序上第 k 与 k+1 个残基之间的真实断口。

    两个判据：残基号跳变、主链 C–N 距离远超 1.33 Å。缺原子时只信编号。
    """
    res: list[tuple[int, dict]] = []
    for ln in lines:
        try:
            num = int(ln[22:26])
        except ValueError:
            continue
        an = ln[12:16].strip()
        if not res or res[-1][0] != num:
            res.append((num, {}))
        if an in ("C", "N"):
            try:
                res[-1][1][an] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError:
                pass

    out = []
    for k in range(len(res) - 1):
        n1, a = res[k]
        n2, b = res[k + 1]
        gap_num = n2 - n1 != 1
        gap_geo = False
        if "C" in a and "N" in b:
            gap_geo = math.dist(a["C"], b["N"]) > PEPTIDE_MAX
        if gap_num or gap_geo:
            out.append(k)
    return out


def _map_by_numbering(seqres, observed):
    """SEQRES 下标 → 残基号：靠编号差值投票，要求一致率够高。"""
    nums = [n for n, _ in observed]
    names = [nm for _, nm in observed]
    pairs = _lcs_pairs(seqres, names)
    if not pairs:
        return None, ["SEQRES 与晶体序列完全没有可比对的部分"]
    off, n_vote = Counter(nums[j] - i for i, j in pairs).most_common(1)[0]
    conf = n_vote / len(pairs)

    mapping, mismatch = {}, []
    for num, nm in observed:
        i = num - off
        if not (0 <= i < len(seqres)):
            continue
        mapping[i] = num
        if not _same_res(seqres[i], nm):
            mismatch.append((num, nm, seqres[i]))
    agree = len(mapping) - len(mismatch)
    if conf < 0.6 or len(mapping) < 0.9 * len(observed) or agree < 0.7 * len(observed):
        return None, [f"编号法不成立（偏移 {off}，投票一致率 {conf:.0%}，"
                      f"残基名相符 {agree}/{len(observed)}）"]
    notes = [f"按编号对位：偏移 {off}，{agree}/{len(observed)} 个残基名相符"]
    return mapping, notes


def _map_by_segments(seqres, observed, breaks):
    """把链按断口切段，每段在 SEQRES 里滑窗找位，段间的间隔就是要补的残基。"""
    segs, s = [], 0
    for k in breaks:
        segs.append((s, k + 1))
        s = k + 1
    segs.append((s, len(observed)))

    mapping, notes, pos = {}, [], 0
    for a, b in segs:
        names = [nm for _, nm in observed[a:b]]
        L = len(names)
        if L == 0:
            continue
        best_ident, best_p = -1, None
        for p in range(pos, len(seqres) - L + 1):
            ident = sum(1 for k in range(L) if _same_res(seqres[p + k], names[k]))
            if ident > best_ident:
                best_ident, best_p = ident, p
        if best_p is None or best_ident < 0.7 * L:
            return None, [f"晶体 {observed[a][0]}–{observed[b - 1][0]} 这 {L} 个残基在 "
                          f"SEQRES 里找不到对应位置（最好只有 {best_ident} 个匹配）"]
        for k in range(L):
            mapping[best_p + k] = observed[a + k][0]
        notes.append(f"按片段对位：晶体 {observed[a][0]}–{observed[b - 1][0]} "
                     f"({L}) → SEQRES {best_p + 1}–{best_p + L}"
                     f"{'' if best_ident == L else f'，{L - best_ident} 个残基名不符'}")
        pos = best_p + L
    return mapping, notes


def build_alignment(seqres, observed, lines, build_termini: bool = False):
    """建 PIR 比对，全部按"位置号"（目标序列第几个残基，1 起）表达。

    返回 (info, error)。info 里有：
        full/tpl   -- 目标序列、模板序列（缺口处 '-'），均单字母
        blocks     -- 要补的缺口 [(首位置, 末位置)]，1 起
        pos2num    -- 位置号 → 原残基号（要新补的为 None）
        first_pos / last_pos -- 模板在比对里的位置范围（structureX 头要用）
        notes      -- 对位方式与可疑之处
    """
    # 断口优先：几何和编号跳变是文件里真实存在的物理断裂，比"编号差值"可靠得多。
    # 编号法只用来交叉验证，两者打架时以断口为准。
    breaks = _break_positions(lines)
    notes = [f"检测到 {len(breaks)} 处主链断口"] if breaks else []
    mapping, notes_seg = _map_by_segments(seqres, observed, breaks)
    if mapping is None:
        m2, notes2 = _map_by_numbering(seqres, observed)
        if m2 is None:
            return None, "；".join(notes_seg + notes2)
        notes.append("几何断口法不成立，退回编号法")
        mapping = m2
    else:
        notes += notes_seg
        m2, _n2 = _map_by_numbering(seqres, observed)
        if m2:
            diff = len(set(mapping) ^ set(m2))
            if diff > 0.05 * len(observed):
                notes.append(f"编号法与断口法有 {diff} 个残基对位不同，以断口法为准")

    idxs = sorted(mapping)
    lo_i = 0 if build_termini else idxs[0]
    hi_i = len(seqres) if build_termini else idxs[-1] + 1
    num2name = {}
    for n, nm in observed:
        num2name.setdefault(n, nm)

    full_chars, tpl_chars, pos2num, unknown = [], [], [], []
    for i in range(lo_i, hi_i):
        aa = _aa1(seqres[i])
        if aa is None:
            unknown.append(seqres[i])
            aa = "X"
        if i in mapping:
            # 模板串写晶体里的真实残基名，写错 MODELLER 会 SequenceMismatchError。
            # 目标串在有坐标的位置照抄模板：SEQRES 和晶体序列不一致的地方（点突变、
            # 缺失的残基数对不上）就不当突变处理，只补缺口，别把已有结构改掉。
            tca = _aa1(num2name[mapping[i]]) or aa
            tpl_chars.append(tca)
            full_chars.append(tca)
            pos2num.append(mapping[i])
        else:
            full_chars.append(aa)
            tpl_chars.append("-")
            pos2num.append(None)
    full, tpl = "".join(full_chars), "".join(tpl_chars)

    if unknown:
        notes.append(f"SEQRES 里有非标准残基 {sorted(set(unknown))}，写作 X 处理")

    # 缺口按位置号给（1 起）；模板串里每个残基名占 3 个字符，正则要扫位置空间
    mask = "".join("." if n is not None else "-" for n in pos2num)
    blocks = [(m.start() + 1, m.end()) for m in re.finditer(r"-+", mask)]

    for a, b in blocks:
        if 1 < a and b < len(pos2num):
            room = pos2num[b] - pos2num[a - 2] - 1      # 原编号留了几个位置
            n = b - a + 1
            if room != n:
                notes.append(f"缺口 位置{a}–{b}（{n} 个残基）：原编号只留了 {room} 个位置")

    if not blocks:
        return None, "没有缺口，不用补"

    first_pos = next(p for p, n in enumerate(pos2num, 1) if n is not None)
    # reversed 里第 r 个是已知残基 → 它的位置号是 len - r
    last_pos = len(pos2num) - next(r for r, n in enumerate(reversed(pos2num))
                                   if n is not None)
    return {"full": full, "tpl": tpl, "blocks": blocks, "pos2num": pos2num,
            "first_pos": first_pos, "last_pos": last_pos, "notes": notes}, None


def _pir(seq: str, width: int = 60) -> str:
    # PIR 规定每条序列必须以 '*' 结尾，少一个 MODELLER 直接报 FileFormatError
    seq = seq.rstrip("*") + "*"
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_alignment(path: Path, code: str, info: dict) -> None:
    body = [
        f">P1;{code}",
        # 模板 PDB 已经按位置号重编，所以这里的范围就是模板首末残基的位置号
        f"structureX:{code}:{info['first_pos']}:A:{info['last_pos']}:A:"
        f"STRUCTURE:ORGANISM:0.00:0.00",
        _pir(info["tpl"]),
        f">P1;{code}_fill",
        # 目标序列一律从 1 编号：MODELLER 不认这里的起始号，写 136 它照样按 1..N 编；
        # 写成 1 反而两边一致，select_atoms 里的号就是位置号
        f"sequence:{code}_fill:1:A:{len(info['full'])}:A::::",
        _pir(info["full"]),
        "",
    ]
    path.write_text("\n".join(body), encoding="utf-8")


def write_script(path: Path, code: str, blocks, n_models: int) -> None:
    sel = ",\n            ".join(
        f"self.residue_range('{a}:A', '{b}:A')" for a, b in blocks)
    path.write_text(f"""from modeller import *
from modeller.automodel import *


class LoopFill(AutoModel):
    def select_atoms(self):
        # 只挑缺口那几段参与建模，模板里已有的残基一个原子都不动
        return Selection(
            {sel},
        )


env = Environ()
env.io.atom_files_directory = ['.']
a = LoopFill(env, alnfile='alignment.ali', knowns='{code}', sequence='{code}_fill')
a.starting_model = 1
a.ending_model = {n_models}
a.make()
""", encoding="utf-8")


def renumber_template(lines: list[str], pos_of_num: dict[int, int]) -> tuple[list[str], int]:
    """模板 PDB 的残基号改成比对位置号，让 structureX 的范围严格对得上。"""
    out, dropped = [], 0
    for ln in lines:
        try:
            num = int(ln[22:26])
        except ValueError:
            out.append(ln)
            continue
        p = pos_of_num.get(num)
        if p is None:
            dropped += 1     # 对不上位的残基不进模板（重号等）
            continue
        out.append(f"{ln[:22]}{p:4d} {ln[27:]}")   # 顺手清掉插入码
    return out, dropped


_MODEL_RE = re.compile(r"(\S+\.pdb)\s+(-?\d+(?:\.\d+)?)")


def pick_best(work: Path, code: str, log: str) -> tuple[Path | None, float, str]:
    """按 molpdf 挑最好的模型：先看日志里的评估表，再退回模型头部的 REMARK。"""
    scores: dict[str, float] = {}
    idx = log.find("Model assessment by MODELLER objective function")
    if idx >= 0:
        for name, val in _MODEL_RE.findall(log[idx:]):
            scores[Path(name).name] = float(val)
    if not scores:
        for f in sorted(work.glob(f"{code}_fill.B*.pdb")):
            m = re.search(r"REMARK\s+6\s+MODELLER OBJECTIVE FUNCTION:\s*(-?[\d.]+)",
                          f.read_text(errors="replace"))
            if m:
                scores[f.name] = float(m.group(1))
    cands = sorted(work.glob(f"{code}_fill.B*.pdb"))
    if not cands:
        return None, 0.0, "没有生成任何模型"
    if not scores:
        return cands[0], 0.0, "日志里没找到 molpdf，取第一个模型"
    best = work / min(scores, key=lambda k: scores[k])
    if not best.exists():
        return cands[0], 0.0, "日志里最好的模型文件不在，取第一个"
    return best, scores[best.name], ""


def _assign_numbers(pos2num: list[int | None]) -> tuple[list[int], bool]:
    """给新补的残基编上号：能塞进原编号的空档就塞，塞不下就顺排并报 False。"""
    out = list(pos2num)
    known = [p for p, n in enumerate(out) if n is not None]
    if not known:
        return [p + 1 for p in range(len(out))], True
    ok, shift = True, 0
    first = known[0]
    last = out[first]
    for k in range(first):
        out[k] = last - (first - k)
    for p in range(first + 1, len(out)):
        n = pos2num[p]
        if n is None:
            last += 1
            out[p] = last
            continue
        num = n + shift
        if num <= last:          # 原编号没给够位置，后面的整体往后顺延
            shift += last + 1 - num
            num = last + 1
            ok = False
        last = num
        out[p] = num
    return out, ok


def _collapse_jumps(atoms: list[str]) -> tuple[list[str], list[str]]:
    """补完之后如果还有"编号跳了但主链连着"的地方，把后面整体顺排掉。

    BCL2 就是这种：晶体编号 39 之后直接跳到 89，其实中间只缺 8 个残基；补完以后
    原子是连着的，但编号留着 42 个空档，下游工具会以为又缺了 42 个残基。
    """
    res: list[tuple[int, dict]] = []
    for ln in atoms:
        try:
            num = int(ln[22:26])
        except ValueError:
            continue
        an = ln[12:16].strip()
        if not res or res[-1][0] != num:
            res.append((num, {}))
        if an in ("C", "N"):
            try:
                res[-1][1][an] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError:
                pass

    need, notes = False, []
    for k in range(len(res) - 1):
        n1, a = res[k]
        n2, b = res[k + 1]
        if n2 - n1 <= 1:
            continue
        d = math.dist(a["C"], b["N"]) if "C" in a and "N" in b else 99.0
        if d <= PEPTIDE_MAX:
            need = True
            notes.append(f"编号 {n1}→{n2} 跳 {n2 - n1 - 1} 个号但主链连着，已顺排")
    if not need:
        return atoms, notes

    order, idx = [], {}
    for ln in atoms:
        try:
            num = int(ln[22:26])
        except ValueError:
            continue
        if num not in idx:
            idx[num] = len(order)
            order.append(num)
    base = order[0]
    out = []
    for ln in atoms:
        try:
            num = int(ln[22:26])
        except ValueError:
            out.append(ln)
            continue
        out.append(f"{ln[:22]}{base + idx[num]:4d} {ln[27:]}")
    return out, notes


def _model_lines(model: Path, chain: str, pos2num: list[int | None]) -> tuple[list[str], list[str]]:
    """读模型，链号改回去，位置号换回原残基号。"""
    atoms = [x for x in model.read_text(errors="replace").splitlines()
             if x.startswith("ATOM")]
    if not atoms:
        return [], []
    try:
        first = int(atoms[0][22:26])
    except ValueError:
        first = 1
    nums, ok = _assign_numbers(pos2num)
    out = []
    for x in atoms:
        try:
            p = int(x[22:26]) - first      # 0 起的位置下标
        except ValueError:
            out.append(x[:21] + chain + x[22:])
            continue
        num = nums[p] if 0 <= p < len(nums) else p + 1
        out.append(f"{x[:21]}{chain}{num:4d} {x[27:]}")
    out, notes = _collapse_jumps(out)
    if not ok:
        notes.insert(0, "补的残基比原编号空档多，跳号之后已整体顺延")
    return out, notes


def graft(src: Path, out: Path, modeled: dict[str, tuple[Path, list]]) -> list[str]:
    """把补好的链拼回原结构：只替换该链的 ATOM 行，其余原样保留。"""
    out_lines, used, notes = [], set(), []
    for ln in src.read_text(errors="replace").splitlines():
        # 修饰残基常以 HETATM 写在链里（TPO308、OAS530）。它们已经并进模板、由模型
        # 输出一份了，原行必须一起去掉，否则同一残基叠两份，mdrun 直接 inf。
        if (ln.startswith("HETATM") and len(ln) > 21 and ln[21] in modeled
                and (ln[17:20].strip().upper() in HEAVY_ATOMS
                     or ln[17:20].strip().upper() in TO_STANDARD)):
            continue
        if ln.startswith("ATOM") and len(ln) > 21 and ln[21] in modeled:
            ch = ln[21]
            if ch in used:
                continue
            used.add(ch)
            mpath, pos2num = modeled[ch]
            body, msgs = _model_lines(mpath, ch, pos2num)
            notes += [f"链 {ch}：{m}" for m in msgs]
            out_lines.extend(body)
            out_lines.append("TER")
            continue
        out_lines.append(ln)
    out.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return notes


def run_chain(py: Path, pdb: Path, chain: str, outdir: Path, args, fasta=None):
    lines, observed, had_icode = read_chain(pdb, chain)
    if not observed:
        _say(f"  链 {chain}：没有蛋白残基，跳过")
        return None
    if had_icode:
        _say(f"  链 {chain}：有插入码（70A 这类），按顺序重编号后再建模")

    seqres = _seqres_by_chain(pdb).get(chain)
    if not seqres and fasta:
        by_chain, default = fasta
        seqres = by_chain.get(chain.upper()) or default
        if seqres:
            _say(f"  链 {chain}：文件里没有 SEQRES，用 --seq 提供的序列（{len(seqres)} 残基）")
    if not seqres:
        _say(f"  链 {chain}：没有 SEQRES，补不了（MODELLER 需要完整序列）。"
             f"加 --seq 全长序列.fasta 即可")
        return None

    info, err = build_alignment(seqres, observed, lines, args.build_termini)
    if err and "没有缺口" in err:
        # 没缺口不代表不用管：CASP3 的编号带 156A 这种插入码，下游按整数编号读
        # 就会以为 157–160 四个残基没了。这种情况把链顺排一遍再原样输出。
        nums = [n for n, _ in observed]
        jumped = any(b - a != 1 for a, b in zip(nums, nums[1:]))
        if not (had_icode or jumped):
            _say(f"  链 {chain}：没有缺口，不用补")
            return "nogap"
        work = outdir / f"chain_{chain}"
        work.mkdir(parents=True, exist_ok=True)
        pos_of_num = {num: p for p, num in enumerate(nums, 1)}
        tpl_lines = [ln[:21] + "A" + ln[22:] for ln in lines]
        tpl_lines, _ = renumber_template(tpl_lines, pos_of_num)
        src = work / "tpl_fill.pdb"
        src.write_text("\n".join(tpl_lines + ["END"]) + "\n")
        _say(f"  链 {chain}：没有缺口，但编号要顺排（{nums[0]}–{nums[0] + len(nums) - 1}）")
        return src, chain, [], nums
    if err:
        _say(f"  链 {chain}：{err}")
        return None

    for n in info["notes"]:
        _say(f"  链 {chain}：· {n}")

    blocks = info["blocks"]
    long_blocks = [(a, b) for a, b in blocks if b - a + 1 > args.max_loop]
    for a, b in blocks:
        tag = "[补]" if b - a + 1 <= args.max_loop else "[补·长]"
        _say(f"  链 {chain}：{tag} 位置 {a}–{b}（{b - a + 1} 个残基）")
    if long_blocks and args.strict:
        _say(f"  链 {chain}：缺口超过 --max-loop {args.max_loop}，"
             "按 --strict 跳过（这种长度补出来构象不可信，建议换更完整的结构）")
        return None
    if long_blocks:
        _say(f"  链 {chain}：! 有缺口超过 {args.max_loop} 个残基，"
             "MODELLER 仍会建，但这些 loop 的构象不可信，用前请看一眼")

    work = outdir / f"chain_{chain}"
    work.mkdir(parents=True, exist_ok=True)

    # 模板只放这条链，链号统一改成 A、残基号改成位置号
    pos_of_num = {}
    for p, num in enumerate(info["pos2num"], 1):
        if num is not None:
            pos_of_num.setdefault(num, p)
    tpl_lines = [ln[:21] + "A" + ln[22:] for ln in lines]
    tpl_lines, dropped = renumber_template(tpl_lines, pos_of_num)
    if dropped:
        _say(f"  链 {chain}：{dropped} 个残基在同一编号下重复，未进模板")
    (work / "tpl.pdb").write_text("\n".join(tpl_lines + ["END"]) + "\n")

    write_alignment(work / "alignment.ali", "tpl", info)
    write_script(work / "run_model.py", "tpl", blocks, args.models)

    _say(f"  链 {chain}：MODELLER 建 {args.models} 个模型（单线程，可能几分钟）…")
    r = subprocess.run([str(py), "run_model.py"], cwd=work, capture_output=True,
                       text=True, encoding="utf-8", errors="replace", timeout=7200)
    log = (r.stdout or "") + (r.stderr or "")
    (work / "modeller.log").write_text(log, encoding="utf-8")
    if r.returncode != 0:
        tail = "\n".join(log.strip().splitlines()[-10:])
        _say(f"  链 {chain}：MODELLER 失败\n{tail}")
        return None

    best, score, note = pick_best(work, "tpl", log)
    if not best:
        _say(f"  链 {chain}：{note}，看 {work / 'modeller.log'}")
        return None
    _say(f"  链 {chain}：选中 {best.name}（molpdf {score:.1f}）"
         f"{'（' + note + '）' if note else ''}")
    return best, chain, long_blocks, info["pos2num"]


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    ap = argparse.ArgumentParser(description="MODELLER 补缺失 loop（晶体结构原位不动）")
    ap.add_argument("pdb")
    ap.add_argument("--out", help="输出 PDB；默认 <输入>_modeled.pdb")
    ap.add_argument("--chain", help="只补这条链，逗号分隔；默认全补")
    ap.add_argument("--models", type=int, default=5, help="每个缺口建几个模型（默认 5）")
    ap.add_argument("--max-loop", type=int, default=30,
                    help="超过这么多残基的缺口给出警告（默认 30）")
    ap.add_argument("--strict", action="store_true", help="有超长缺口就整条链跳过")
    ap.add_argument("--build-termini", action="store_true", help="首尾悬空段也建")
    ap.add_argument("--seq", action="append", metavar="FASTA",
                    help="完整序列（FASTA）。文件里没有 SEQRES 时用这个补——PyMOL/AutoDock "
                         "清洗过的受体常没有 SEQRES。记录名单字母（A/B）按链分配，"
                         "否则当作默认序列给所有链用")
    ap.add_argument("--keep-work", action="store_true", help="保留中间文件")
    a = ap.parse_args(argv)

    if os.name == "nt":
        try:
            import ctypes
            ctypes.windll.kernel32.SetPriorityClass(
                ctypes.windll.kernel32.GetCurrentProcess(), 0x00004000)
        except Exception:
            pass

    py = find_modeller_python()
    if not py:
        sys.exit("没找到装 MODELLER 的环境。先建：\n"
                 "  conda create -n modeller -c salilab -c conda-forge python=3.12 modeller -y\n"
                 "再在 envs/modeller/Library/modeller/modlib/modeller/config.py 里填 license")

    src = Path(a.pdb)
    if not src.exists():
        sys.exit(f"找不到 {src}")
    out = Path(a.out) if a.out else src.with_name(src.stem + "_modeled.pdb")

    chains = ([c.strip() for c in a.chain.split(",")] if a.chain
              else sorted({ln[21] for ln in src.read_text(errors='replace').splitlines()
                           if ln.startswith("ATOM")}))
    _say(f"处理 {src.name}：链 {chains}")
    fasta = _load_fasta(a.seq) if a.seq else None

    workroot = Path(tempfile.mkdtemp(prefix="modeller_", dir=str(out.parent)))
    try:
        modeled: dict[str, tuple[Path, list]] = {}
        long_warn: dict[str, list] = {}
        nogap = 0
        for ch in chains:
            res = run_chain(py, src, ch, workroot, a, fasta)
            if res == "nogap":
                nogap += 1
            elif res:
                best, ch2, longs, pos2num = res
                modeled[ch2] = (best, pos2num)
                if longs:
                    long_warn[ch2] = longs

        if not modeled:
            if nogap == len(chains):
                _say("\n所有链都没有缺口，不需要补。")
                return 0
            _say("\n没有任何链补成功。")
            return 1

        notes = graft(src, out, modeled)
        _say(f"\n写出 {out}")
        for n in notes:
            _say(f"! {n}")

        _say("复查：")
        try:
            import check_protein
            rep = check_protein.analyse(out)
            if rep["issues"]:
                _say(f"  体检 {rep['verdict']}：" + "；".join(
                    i["what"] for i in rep["issues"][:4]))
            else:
                _say("  体检 PASS：结构完整")
        except Exception as e:
            _say(f"  （体检没跑成：{e}）")
            _say(f'  手工复查：python scripts/check_protein.py "{out}"')
        if long_warn:
            _say("! 以下长 loop 是 MODELLER 硬补的，构象可信度低：")
            for ch, bl in long_warn.items():
                _say(f"  链 {ch}：" + "、".join(f"位置 {x}–{y}" for x, y in bl))
        return 0
    finally:
        if not a.keep_work:
            shutil.rmtree(workroot, ignore_errors=True)
        else:
            _say(f"中间文件：{workroot}")


if __name__ == "__main__":
    sys.exit(main())
