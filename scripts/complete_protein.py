#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""补全残缺蛋白：建模缺失的 loop、补齐缺失侧链重原子、替换非标准残基。

体检工具 check_protein.py 报 FAIL 的，通常就是这个脚本的活。用的是 PDBFixer
（OpenMM 生态，BSD 许可，无需注册）——它从模板库取残基、建 loop，再可选地在
隐式溶剂里做一次带约束的能量最小化把 clash 消掉。

注意：补出来的 loop 是"一个合理构象"，不是唯一解。无序区越长越不可信，所以默认
只补 ≤ 15 个残基的内部缺口，首尾悬空段一律不补（那是本来就无序的尾巴，补了只会
在 MD 里乱飘）。要改就调 --max-loop / --build-termini。

**只用 CPU，绝不碰 GPU。** OpenMM 默认会挑最快的平台（本机是 CUDA），那样最小化
会跟正在跑的 mdrun 抢显卡、把主线吞吐拉下来。所以这里强制 CPU 平台，并且在进程里
把 CUDA 设备整个屏蔽掉；算一个 500 残基的蛋白只要几十秒，不值得动 GPU。

装环境（一次性，约 5 分钟）：
    conda create -n mdprep -c conda-forge python=3.11 openmm pdbfixer -y
    C:\\Users\\<你>\\miniconda3\\envs\\mdprep\\python.exe scripts\\complete_protein.py ...

用法：
    python scripts/complete_protein.py rec.pdb --out rec_fixed.pdb --dry-run
    python scripts/complete_protein.py rec.pdb --out rec_fixed.pdb
    python scripts/complete_protein.py rec.pdb --out rec_fixed.pdb --chains A --max-loop 8
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

H_Y_NAMES = ("H", "D")


def _say(msg=""):
    print(msg, flush=True)


def load_fixer(path: Path):
    from pdbfixer import PDBFixer
    return PDBFixer(filename=str(path))


def _seqres_by_chain(path: Path) -> dict[str, list[str]]:
    """从 PDB 里读 SEQRES，按链给出完整序列。"""
    out: dict[str, list[str]] = {}
    last_serial: dict[str, int] = {}
    for line in Path(path).read_text(errors="replace").splitlines():
        if not line.startswith("SEQRES"):
            continue
        ch = line[11] if len(line) > 11 else " "
        try:
            serial = int(line[7:10])
        except ValueError:
            serial = None
        # 有些文件里 SEQRES 整块重复了几遍（STAT3.pdb 就是 3 遍），
        # 序号一往回走就说明换了一遍，后面的不要
        if serial is not None and ch in last_serial and serial <= last_serial[ch]:
            continue
        if serial is not None:
            last_serial[ch] = serial
        for i in range(19, 70, 4):
            nm = line[i:i + 3].strip()
            if nm:
                out.setdefault(ch, []).append(nm)
    return out


def _lcs_pairs(a: list[str], b: list[str]) -> list[tuple[int, int]]:
    """最长公共子序列，返回配对的 (i, j)。用来把 SEQRES 对到实际残基上。"""
    n, m = len(a), len(b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n - 1, -1, -1):
        ai = a[i]
        row, nxt = dp[i], dp[i + 1]
        for j in range(m - 1, -1, -1):
            row[j] = nxt[j + 1] + 1 if ai == b[j] else (nxt[j] if nxt[j] >= row[j + 1] else row[j + 1])
    pairs, i, j = [], 0, 0
    while i < n and j < m:
        if a[i] == b[j]:
            pairs.append((i, j)); i += 1; j += 1
        elif dp[i + 1][j] >= dp[i][j + 1]:
            i += 1
        else:
            j += 1
    return pairs


# SEQRES 里偶尔会出现修饰残基（OAS/SEP/TPO/MSE…），PDBFixer 没有对应的模板，
# 直接建会崩。按最常见的酸/酰胺母体替换——这和 Amber/CHARMM 处理修饰残基的
# 惯例一致：保留骨架，丢掉修饰基团。
MODIFIED_TO_STD = {
    "MSE": "MET", "FME": "MET", "SEP": "SER", "TPO": "THR", "PTR": "TYR",
    "CSO": "CYS", "CSD": "CYS", "CME": "CYS", "OCY": "CYS", "CSS": "CYS",
    "KCX": "LYS", "LLP": "LYS", "MLY": "LYS", "M3L": "LYS", "ALY": "LYS",
    "HYP": "PRO", "PCA": "GLU", "OAS": "SER", "OAT": "ALA", "NLE": "LEU",
    "SAR": "GLY", "ABA": "ALA", "DAL": "ALA", "DLE": "LEU", "DVA": "VAL",
}


def _derive_missing(fixer, seqres: dict[str, list[str]]) -> dict:
    """自己做序列比对，算出 {（链下标, 插入位置）: [残基名…]}。

    PDBFixer 自带的 findMissingResidues 要求 SEQRES 和实际残基严格按序对齐，
    一旦编号重复（同一编号出现两次）或整体错位就直接放弃、一个缺口都不报。
    CASP3 / PTGS2 就是这种情况：明明有 SEQRES，它却说"没有编号缺口"。
    """
    out = {}
    for chain in fixer.topology.chains():
        names = seqres.get(chain.id)
        if not names:
            continue
        obs = [r.name for r in chain.residues()]
        if not obs:
            continue
        pairs = _lcs_pairs(names, obs)
        if len(pairs) < 0.9 * len(obs):        # 对不上的链不要硬凑
            continue
        prev_i = -1
        for (i, j) in pairs + [(len(names), len(obs))]:
            if i - prev_i > 1:                  # 中间这一段在 SEQRES 里有、结构里没有
                block = [MODIFIED_TO_STD.get(nm, nm) for nm in names[prev_i + 1:i]]
                bad = [nm for nm in block if nm not in _STD_AA]
                if bad:
                    print(f'  ! 链 {chain.id} 这一段缺口里有认不出的残基 {set(bad)}，'
                          f'跳过不补', flush=True)
                else:
                    out[(chain.index, j)] = block
            prev_i = i
    return out


_STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def describe_gaps(fixer):
    """把 PDBFixer 找到的缺口整理成人能读的列表。"""
    topo = fixer.topology
    chains = list(topo.chains())
    rows = []
    for (ci, ri), resnames in sorted(fixer.missingResidues.items()):
        chain = chains[ci]
        reslist = list(chain.residues())
        before = reslist[ri - 1] if ri > 0 else None
        after = reslist[ri] if ri < len(reslist) else None
        where = f'{before.name}{before.id} → {after.name}{after.id}' if (before and after) else (
            f'N 端延伸到 {after.name}{after.id}' if after else f'C 端延伸到 {before.name}{before.id}')
        rows.append({"chain": chain.id, "key": (ci, ri), "n": len(resnames),
                     "where": where, "resnames": resnames,
                     "internal": bool(before and after)})
    return rows


def main(argv=None):
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    ap = argparse.ArgumentParser(description="用 PDBFixer 补全缺失 loop / 缺失原子")
    ap.add_argument("pdb", help="输入 PDB")
    ap.add_argument("--out", help="输出 PDB（不给就和输入同名加 _fixed）")
    ap.add_argument("--chains", help="只保留这些链，逗号分隔（如 A 或 A,B）；默认全留")
    ap.add_argument("--seqres-from", metavar="原始PDB",
                    help="从这个文件借 SEQRES 记录。PyMOL/受体处理工具清洗过的 *_clean.pdb "
                         "常常没有 SEQRES，PDBFixer 就看不到缺口；借原版 PDB 的 SEQRES 即可")
    ap.add_argument("--max-loop", type=int, default=15,
                    help="内部缺口最多补多少个残基，超过就跳过（默认 15）")
    ap.add_argument("--build-termini", action="store_true",
                    help="首尾悬空段也补（默认不补）")
    ap.add_argument("--keep-het", action="store_true",
                    help="保留非蛋白组分（配体/离子/水）；默认全删")
    ap.add_argument("--keep-water", action="store_true", help="配合 --keep-het 保留水")
    ap.add_argument("--no-minimize", action="store_true",
                    help="不做约束最小化（补出来的 loop 可能和邻位有 clash）")
    ap.add_argument("--min-steps", type=int, default=300, help="最小化步数（默认 300）")
    ap.add_argument("--threads", type=int, default=1,
                    help="最小化用几个 CPU 线程（默认 1，避免抢 MD 的核）")
    ap.add_argument("--use-gpu", action="store_true",
                    help="允许用 GPU 做最小化（默认禁止：会和 mdrun 抢显卡，别在跑 MD 时开）")
    ap.add_argument("--normal-priority", action="store_true",
                    help="用正常进程优先级（默认降到 BelowNormal，让 mdrun 优先）")
    ap.add_argument("--dry-run", action="store_true", help="只报告要补什么，不写文件")
    a = ap.parse_args(argv)

    # 必须在 import openmm 之前屏蔽 CUDA，否则 OpenMM 会自己挑到 GPU 平台上，
    # 和正在跑的 mdrun 抢显卡
    import os
    if not a.use_gpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    # 降到 BelowNormal：CPU 抢的时候永远让 mdrun 先跑
    if os.name == "nt" and not a.normal_priority:
        try:
            import ctypes
            BELOW_NORMAL_PRIORITY_CLASS = 0x00004000
            ctypes.windll.kernel32.SetPriorityClass(
                ctypes.windll.kernel32.GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS)
        except Exception:
            pass

    try:
        import openmm  # noqa: F401
    except ImportError:
        sys.exit("没装 openmm/pdbfixer。先建环境：\n"
                 "  conda create -n mdprep -c conda-forge python=3.11 openmm pdbfixer -y\n"
                 "然后用 miniconda3/envs/mdprep/python.exe 跑这个脚本。")

    src = Path(a.pdb)
    if not src.exists():
        sys.exit(f"找不到文件：{src}")
    out = Path(a.out) if a.out else src.with_name(src.stem + "_fixed.pdb")

    if a.seqres_from:
        src = _inject_seqres(src, Path(a.seqres_from))

    _say(f"读取 {src.name} …")
    fixer = load_fixer(src)
    topo = fixer.topology

    n_before = sum(1 for r in topo.residues() if r.name not in ("HOH", "WAT"))
    _say(f"  原始：{len(list(topo.chains()))} 条链，{n_before} 个残基")

    keep_chains = None
    if a.chains:
        keep_chains = {c.strip() for c in a.chains.split(",") if c.strip()}

    # ---- 非标准残基（MSE 等）换成标准残基
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        _say("  替换非标准残基：" + "、".join(r.name for r, _ in fixer.nonstandardResidues))
        fixer.replaceNonstandardResidues()

    # ---- 杂原子
    if not a.keep_het:
        het = {}
        for r in fixer.topology.residues():
            if r.name not in _POLY and r.name not in ("HOH", "WAT"):
                het[r.name] = het.get(r.name, 0) + 1
        if het:
            _say("  删除杂原子：" + "、".join(f"{k}×{v}" for k, v in het.items()))
        fixer.removeHeterogens(keepWater=a.keep_water)

    # ---- 只留指定链
    if keep_chains:
        drop = [c for c in list(fixer.topology.chains()) if c.id not in keep_chains]
        if drop:
            fixer.removeChains([c.index for c in drop])
            _say(f"  保留链 {sorted(keep_chains)}，删掉 {len(drop)} 条")

    # ---- 找缺口。必须放在删杂原子、删链之后：缺口的 key 里带链下标，
    # 拓扑一变下标就错位，补出来的残基会插到别的链上去。
    fixer.findMissingResidues()
    seqres = _seqres_by_chain(src)
    derived = _derive_missing(fixer, seqres) if seqres else {}
    if not fixer.missingResidues and derived:
        fixer.missingResidues = derived
        _say(f"  PDBFixer 的 SEQRES 对齐失败，改用序列比对补出 {len(derived)} 处缺口")
    gaps = describe_gaps(fixer)

    chosen, skipped = [], []
    for g in gaps:
        why = None
        if not g["internal"] and not a.build_termini:
            why = "首尾悬空段（默认不补）"
        elif g["n"] > a.max_loop:
            why = f'长度 {g["n"]} > --max-loop {a.max_loop}'
        (skipped if why else chosen).append((g, why))

    if gaps:
        _say(f"  发现 {len(gaps)} 处缺口：补 {len(chosen)} 处，跳过 {len(skipped)} 处")
        for g, why in chosen:
            _say(f'    [补] 链 {g["chain"]}  {g["where"]}  {g["n"]} 个残基')
        for g, why in skipped:
            _say(f'    [跳] 链 {g["chain"]}  {g["where"]}  {g["n"]} 个残基 —— {why}')
    else:
        _say("  没有编号缺口")
        # 这时体检报告里的缺口就补不了，必须说清楚，别让人以为已经修好了
        try:
            from check_protein import analyse
            rep = analyse(src)
            missed = sum(len(c["gaps_numbering"]) for c in rep["chains"])
            if missed:
                _say(f"  ! 体检发现 {missed} 处编号缺口，但推不出缺的到底是哪几个残基"
                     "（没有 SEQRES 或序列对不上）。")
                _say("    这些位置补不了，需要手工处理或换 MODELLER / SWISS-MODEL 建模。")
        except Exception:
            pass

    total_new = sum(g["n"] for g, _ in chosen)
    if a.dry_run:
        _say(f"\n[dry-run] 会补 {total_new} 个残基，未写文件。")
        return 0

    # ---- 按选择裁剪要建的缺口
    fixer.missingResidues = {g["key"]: g["resnames"] for g, _ in chosen}

# ---- 记下"原有残基"。必须用 (链ID, 残基ID) 做键，不能用 index：
    # addMissingResidues 会往链中间插残基，插完之后所有残基的 index 都往后挪，
    # 拿旧 index 去加位置约束会把约束加到错的原子身上，把结构整个拉变形。
    original_keys = {(c.id, r.id)
                     for c in fixer.topology.chains() for r in c.residues()}

    # ---- 补残基、补原子
    if total_new:
        _say(f"  建模 {total_new} 个缺失残基 …")
    fixer.findMissingAtoms()
    n_missing_atoms = sum(len(v) for v in fixer.missingAtoms.values())
    if n_missing_atoms:
        _say(f"  补齐 {n_missing_atoms} 个缺失原子（{len(fixer.missingAtoms)} 个残基）")
    fixer.addMissingAtoms()

    _say("  加氢（pH 7.0）…")
    fixer.addMissingHydrogens(7.0)

    # ---- 约束最小化，只放开新补的部分
    if not a.no_minimize and (total_new or n_missing_atoms):
        try:
            _say(f"  约束最小化 {a.min_steps} 步（CPU {a.threads} 线程）…")
            # 用 id 重新映射成当前拓扑里的 index
            keep_res = {(c.index, r.index)
                        for c in fixer.topology.chains() for r in c.residues()
                        if (c.id, r.id) in original_keys}
            fixer = relax(fixer, keep_res, steps=a.min_steps, threads=a.threads)
            _say("  最小化完成")
        except Exception as e:                     # 最小化失败不影响输出
            _say(f"  ! 最小化失败（{type(e).__name__}: {e}），输出未减速的结构")

    # ---- 写文件：去掉氢，交给 pdb2gmx 重新加
    from openmm.app import PDBFile
    _say(f"  写出 {out.name} …")
    keep_h = False
    with out.open("w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    _strip_hydrogens(out, keep_h)

    _say(f"完成：{out}")

    # 写完立刻复查，别让人拿着一个还有缺口的结构去跑 MD
    try:
        from check_protein import analyse
        rep = analyse(out)
        errs = [i for i in rep["issues"] if i["level"] == "ERROR"]
        if errs:
            _say(f"\n复查：还有 {len(errs)} 处问题")
            for i in errs[:10]:
                _say(f'  [{i["chain"]}] {i["what"]}')
            if len(errs) > 10:
                _say(f'  …另有 {len(errs) - 10} 处')
            _say("PDBFixer 只能从 SEQRES 建缺口。剩下的这些要手工补，或者换")
            _say("MODELLER / SWISS-MODEL 建模后再回来跑流水线。")
            return 1
        _say("复查：通过，可以进流水线")
    except Exception as e:
        _say(f"（复查跳过：{e}）")
        _say(f'建议自己跑一遍：python scripts/check_protein.py "{out}"')
    return 0


_POLY = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN", "MSE", "SEC", "UNK",
}


def _inject_seqres(src: Path, donor: Path) -> Path:
    """把 donor 的 SEQRES 记录塞进 src，写到临时文件。

    PyMOL / AutoDock 的受体处理工具清洗 PDB 时会把 SEQRES 一起丢掉，PDBFixer
    于是完全看不到中间缺了哪几段、一个 loop 都不建。原版 PDB 里通常还留着。
    """
    all_seqres = [ln for ln in donor.read_text(errors="replace").splitlines()
                  if ln.startswith("SEQRES")]
    if not all_seqres:
        sys.exit(f"{donor.name} 里没有 SEQRES 记录，借不了")
    body = [ln for ln in src.read_text(errors="replace").splitlines() if not ln.startswith("SEQRES")]

    # 只保留目标文件里真有的链。原版 PDB 常有 A/B 多条链，清洗后只剩一条，
    # SEQRES 全带过来会对不上号，PDBFixer 一样建不出 loop。
    chains_in_src = {ln[21] for ln in body if ln.startswith(("ATOM", "HETATM"))}
    seqres = [ln for ln in all_seqres if ln[11] in chains_in_src]
    if not seqres:
        sys.exit(f"{donor.name} 的 SEQRES 链 {sorted({ln[11] for ln in all_seqres})} "
                 f"和目标文件 {sorted(chains_in_src)} 对不上")
    # 每条链的残基数按实际保留的名字重算（列 14-17）
    counts: dict[str, int] = {}
    for ln in seqres:
        ch = ln[11]
        counts[ch] = counts.get(ch, 0) + sum(
            1 for i in range(19, 70, 4) if ln[i:i + 3].strip())
    seqres = [ln[:13] + f"{counts[ln[11]]:>4}" + ln[17:] for ln in seqres]
    idx = next((i for i, ln in enumerate(body)
                if ln.startswith(("ATOM", "HETATM"))), len(body))
    dst = Path(tempfile.mkdtemp()) / (src.stem + "_sq.pdb")
    dst.write_text("\n".join(body[:idx] + seqres + body[idx:]) + "\n")
    return dst


def _strip_hydrogens(path: Path, keep_h: bool) -> None:
    """PDBFixer 写出来的带氢；pdb2gmx 会自己加，去掉免得重复。"""
    if keep_h:
        return
    lines = []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            elem = line[76:78].strip() if len(line) > 76 else ""
            if elem == "H" or name.startswith("H") or (name[:1].isdigit() and name[1:2] == "H"):
                continue
        lines.append(line)
    path.write_text("\n".join(lines) + "\n")


def relax(fixer, keep_res, steps=300, threads=1, k=10.0):
    """真空里带位置约束最小化；原有残基被拉住，新补的 loop 自由落地。

    为什么用真空而不是隐式溶剂：OpenMM 的 CPU 平台跑 GBSA 是 O(N²) 的重活，
    3700 原子的体系一步要 10–30 秒，200 步得跑一小时，完全没法在 MD 旁边并行。
    真空 + 1.2 nm 截断一步 0.05 秒，快 200 倍以上。

    这样做不损失什么：这一步的目的只有一个——把 PDBFixer 从模板拼出来的 loop
    里的原子撞车推开。真正的结构弛豫交给 GROMACS 的 EM/NVT/NPT（显式水、正确的
    力场条件）。原有原子的位置约束（10 kcal/mol/Å²）保证实验结构不会被真空效应
    带偏，只有新补的残基在动。
    """
    import openmm
    from openmm import unit
    from openmm.app import ForceField, Modeller, CutoffNonPeriodic, HBonds

    ff = ForceField("amber14-all.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    system = ff.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic,
                             nonbondedCutoff=1.2 * unit.nanometer,
                             constraints=HBonds, rigidWater=True)

    # 给原有原子加位置约束；新补的残基不加，让它们松弛
    rest = openmm.CustomExternalForce(
        "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    rest.addGlobalParameter("k", k * unit.kilocalorie_per_mole / unit.angstrom**2)
    for p in ("x0", "y0", "z0"):
        rest.addPerParticleParameter(p)

    pos = modeller.positions
    topo = modeller.topology
    for chain in topo.chains():
        for res in chain.residues():
            if (chain.index, res.index) not in keep_res:
                continue
            for atom in res.atoms():
                # x0 必须和 OpenMM 内部单位一致（nm）。以前在这里乘了 10 换成 Å，
                # 于是每个受约束的原子都被往 10 倍远的地方拽，PTGS2 整条链的肽键
                # 被拉长到 2.3 Å，看着像"到处断链"
                p = pos[atom.index].value_in_unit(unit.nanometer)
                rest.addParticle(atom.index, [v for v in p])
    system.addForce(rest)

    # 注：试过在这里补建接缝肽键、也试过放开接缝两侧的约束，两种做法都会让
    # 最小化发散（STAT3 的接缝从 4 Å 恶化到 11 Å），所以都不做。
    # PDBFixer 补出来的接缝键长偏长（2–5 Å），交给 GROMACS 的 EM 收——
    # pdb2gmx 是按残基序列建键的，不看坐标，拓扑本身是连通的。

    integ = openmm.LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                      2 * unit.femtosecond)
    # 显式指定 CPU 平台。不写这一句，OpenMM 会按"最快"去选，本机上就是 CUDA，
    # 于是最小化跑去和 mdrun 抢同一块显卡。
    plat = openmm.Platform.getPlatformByName("CPU")
    print(f"  OpenMM 平台：{plat.getName()}", flush=True)
    ctx = openmm.Context(system, integ, plat, {"Threads": str(max(1, threads))})
    ctx.setPositions(modeller.positions)
    openmm.LocalEnergyMinimizer.minimize(ctx, 10.0, steps)

    newpos = ctx.getState(getPositions=True).getPositions()
    del ctx, integ
    return _Refixed(fixer.topology, newpos)


class _Refixed:
    """最小化只换坐标，拓扑不变。"""
    def __init__(self, topology, positions):
        self.topology = topology
        self.positions = positions


if __name__ == "__main__":
    sys.exit(main())
