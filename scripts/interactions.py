#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""interactions.py —— 配体-蛋白相互作用：接触频率、氢键配对、配体聚类、2D 示意图。

不做逐帧的 gmx 命令堆叠，而是先用 `gmx traj -ox` 把 Protein_LIG 的坐标
一次性导出成文本，再用 numpy 统一算——这样接触、氢键、盐桥都能拿到
“逐对占有率”，gmx hbond 只给总数，给不了这个。

判据（文献里最常用的写法）：
  接触    任意重原子间距 ≤ 4.5 Å
  氢键    D···A ≤ 3.5 Å 且 D–H···A ≥ 120°
  盐桥    带正电残基 N 与带负电残基 O 间距 ≤ 4.0 Å

聚类用 `gmx cluster`（GROMOS 法，cutoff 0.2 nm），以主链叠合后看配体位姿，
能把“结合态/解离态”分开。

2D 示意图用 RDKit 读 `01_ligand/pose_h.mol2` 的键连接，画成 LigPlot 风格的
径向图：配体在中间，相互作用的残基围在四周，蓝虚线氢键、灰点线疏水接触、
红点划线盐桥，括号里是该帧组的占有率。

用法：
    python scripts/interactions.py --job Glc_FABP5 --replicates 1 2
    python scripts/interactions.py --job Glc_FABP5 --cutoff 0.4 --dt 0.5
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))


# --------------------------------------------------------------------------
# 输入准备
# --------------------------------------------------------------------------
def read_gro(path: Path) -> list[tuple[int, str, str]]:
    """读 .gro 前三列：残基号、残基名、原子名。"""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    n = int(lines[1].split()[0])
    atoms = []
    for ln in lines[2:2 + n]:
        atoms.append((int(ln[0:5]), ln[5:10].strip(), ln[10:15].strip()))
    return atoms


def extract_xyz(cfg, m: dict, dt_ns: float = 0.2) -> tuple[np.ndarray, np.ndarray]:
    """返回 (时间 ns, 坐标数组 [帧, 原子, 3]，单位 nm)。"""
    from viz import _gmx
    out = m["out"] / f"coords{m['sfx']}.xvg"
    if not out.exists():
        _gmx(cfg, ["traj", "-s", m["ref"], "-f", m["prot"], "-n", m["ndx"],
                   "-ox", out, "-dt", dt_ns, "-tu", "ns"],
             cwd=m["job"], stdin_text="Protein_LIG\n",
             log=m["out"] / f"_traj{m['sfx']}.log")
    arr = np.loadtxt(out, comments=["#", "@"])
    if arr.ndim == 1:
        arr = arr[None, :]
    t = arr[:, 0]
    xyz = arr[:, 1:].reshape(len(arr), -1, 3)
    return t, xyz


# --------------------------------------------------------------------------
# 相互作用统计
# --------------------------------------------------------------------------
def _element(name: str) -> str:
    """原子名取首字母作元素（OG1→O、NZ→N、HG1→H），够用于蛋白/糖类。"""
    for ch in name:
        if ch.isalpha():
            return ch.upper()
    return "?"


def topology_masks(atoms, resname="LIG"):
    """返回配体重原子、蛋白重原子、供体、受体等下标的索引。"""
    is_lig = np.array([a[1] == resname for a in atoms])
    is_h = np.array([a[2].upper().startswith("H") for a in atoms])
    elem = np.array([_element(a[2]) for a in atoms])
    heavy = ~is_h
    return {
        "is_lig": is_lig, "heavy": heavy, "elem": elem,
        "lig_heavy": np.where(is_lig & heavy)[0],
        "prot_heavy": np.where(~is_lig & heavy)[0],
        "lig_all": np.where(is_lig)[0],
        "prot_all": np.where(~is_lig)[0],
        "lig_o": np.where(is_lig & heavy & np.isin(elem, ["O", "N"]))[0],
        "prot_o": np.where(~is_lig & heavy & np.isin(elem, ["O", "N"]))[0],
    }


def _donors(xyz, atoms, mask, which: str):
    """按最近重原子判断每个 H 的母原子；只保留母原子为 O/N 的供体。"""
    idx = mask[which]
    is_h = np.array([a[2].upper().startswith("H") for a in atoms])
    h_idx = idx[is_h[idx]]
    heavy_idx = idx[~is_h[idx]]
    if len(h_idx) == 0:
        return []
    x0 = xyz[0]
    out = []
    for h in h_idx:
        d = np.linalg.norm(x0[heavy_idx] - x0[h], axis=1)
        j = int(np.argmin(d))
        if d[j] <= 0.13 and _element(atoms[heavy_idx[j]][2]) in ("O", "N"):
            out.append((int(h), int(heavy_idx[j])))
    return out


def hbond_table(xyz, atoms, masks, resname="LIG",
                dmax=0.35, amin=120.0, bound=None):
    """蛋白-配体氢键的逐对占有率。

    `bound` 是逐帧的布尔掩码（配体与蛋白有接触）；给了就额外算
    “只在结合帧里”的占有率——弱结合体系里这个数比全程占有率好读得多。
    """
    donors_p = _donors(xyz, atoms, masks, "prot_all")
    donors_l = _donors(xyz, atoms, masks, "lig_all")
    acc_p, acc_l = masks["prot_o"], masks["lig_o"]
    nf = len(xyz)
    rows: dict[tuple, dict] = {}

    def scan(donors, acceptors, donor_is_lig):
        if not donors or len(acceptors) == 0:
            return
        hi = np.array([d[0] for d in donors])
        di = np.array([d[1] for d in donors])
        for f in range(nf):
            X = xyz[f]
            da = np.linalg.norm(X[di][:, None, :] - X[acceptors][None, :, :], axis=2)
            if not (da <= dmax).any():
                continue
            v1 = X[di][:, None, :] - X[hi][:, None, :]
            v2 = X[acceptors][None, :, :] - X[hi][:, None, :]
            c = (v1 * v2).sum(2) / (np.linalg.norm(v1, axis=2) *
                                    np.linalg.norm(v2, axis=2) + 1e-9)
            ang = np.degrees(np.arccos(np.clip(c, -1, 1)))
            hit = (da <= dmax) & (ang >= amin)
            for a, b in np.argwhere(hit):
                if donor_is_lig:
                    key = (int(acceptors[b]), int(di[a]))
                else:
                    key = (int(di[a]), int(acceptors[b]))
                r = rows.setdefault(key, {"n": 0, "d": [], "ang": []})
                r["n"] += 1
                r["d"].append(float(da[a, b]))
                r["ang"].append(float(ang[a, b]))

    scan(donors_p, acc_l, False)   # 蛋白供体 → 配体受体
    scan(donors_l, acc_p, True)    # 配体供体 → 蛋白受体

    nb = int(np.sum(bound)) if bound is not None else nf
    out = []
    for (p, l), r in rows.items():
        pr, pa, ln, la = atoms[p][0], atoms[p][2], atoms[l][0], atoms[l][2]
        out.append({
            "resnr": pr, "resname": atoms[p][1], "prot_atom": pa,
            "lig_atom": la,
            "occ": r["n"] / nf,
            "occ_bound": r["n"] / max(1, nb),
            "dist": float(np.mean(r["d"])),
            "angle": float(np.mean(r["ang"])),
        })
    out.sort(key=lambda r: -r["occ_bound"])
    return out


def contact_table(xyz, atoms, masks, cutoff=0.45, lig="LIG"):
    """逐残基接触占有率 + 平均最小距离。"""
    from scipy.spatial.distance import cdist
    ph, lh = masks["prot_heavy"], masks["lig_heavy"]
    resnr = np.array([atoms[i][0] for i in range(len(atoms))])[ph]
    resname = np.array([atoms[i][1] for i in range(len(atoms))])[ph]
    atomname = np.array([atoms[i][2] for i in range(len(atoms))])[ph]
    keys = sorted({(int(r), n) for r, n in zip(resnr, resname)})
    lig_names = np.array([atoms[i][2] for i in masks["lig_heavy"]])
    lig_c = np.where(np.isin(masks["elem"][lh], ["C"]))[0]  # 疏水接触挂到碳上
    groups = {k: np.where((resnr == k[0]) & (resname == k[1]))[0] for k in keys}
    hits = {k: 0 for k in keys}
    dmin = {k: [] for k in keys}
    best = {k: (9.9, "") for k in keys}   # 该残基最近的一次接触对应哪个配体原子
    best_c = {k: (9.9, "") for k in keys}
    bound = []
    for f in range(len(xyz)):
        d = cdist(xyz[f][ph], xyz[f][lh])
        amin = d.min(1)
        bound.append(bool(amin.min() <= cutoff))
        for k, idx in groups.items():
            sub = d[idx]
            i, j = np.unravel_index(int(np.argmin(sub)), sub.shape)
            dm = float(sub[i, j])
            dmin[k].append(dm)
            if dm < best[k][0]:
                best[k] = (dm, str(lig_names[j]))
            if len(lig_c):
                subc = sub[:, lig_c]
                ic, jc = np.unravel_index(int(np.argmin(subc)), subc.shape)
                dc = float(subc[ic, jc])
                if dc < best_c[k][0]:
                    best_c[k] = (dc, str(lig_names[lig_c[jc]]))
            if dm <= cutoff:
                hits[k] += 1
    nf = len(xyz)
    nb = max(1, int(np.sum(bound)))
    out = []
    for k in keys:
        idx = groups[k]
        out.append({
            "resnr": k[0], "resname": k[1],
            "occ": hits[k] / nf,
            "occ_bound": hits[k] / nb,
            "min_dist": float(np.min(dmin[k])),
            "mean_dist": float(np.mean(dmin[k])),
            "lig_atom": best[k][1],
            "lig_c_atom": best_c[k][1],
            "atoms": [str(atomname[i]) for i in idx][:6],
        })
    out.sort(key=lambda r: -r["occ_bound"])
    return out, np.array(bound)


def salt_bridges(xyz, atoms, masks, resname="LIG", cutoff=0.40, bound=None,
                 neg_names=None):
    """蛋白正电残基与配体羧基 O 之间的盐桥。

    注意只认真正的羧酸根：葡萄糖的羟基 O 是中性，不该出现在盐桥里。
    `neg_names` 由 `ligand_anionic_atoms()` 从 mol2 连接性判定。
    """
    if not neg_names:
        return []
    lig_neg = [i for i in masks["lig_o"] if atoms[i][2] in neg_names]
    if not lig_neg:
        return []
    pos_sets = {"ARG": ["NH1", "NH2", "NE"], "LYS": ["NZ"],
                "HIS": ["ND1", "NE2"]}
    prot_pos = [i for i in masks["prot_heavy"]
                if atoms[i][1] in pos_sets and atoms[i][2] in pos_sets[atoms[i][1]]]
    if not prot_pos:
        return []
    nf = len(xyz)
    hits = {}
    for f in range(nf):
        X = xyz[f]
        d = np.linalg.norm(X[prot_pos][:, None, :] - X[lig_neg][None, :, :], axis=2)
        for a, b in np.argwhere(d <= cutoff):
            key = (prot_pos[a], lig_neg[b])
            hits[key] = hits.get(key, 0) + 1
    nb = int(np.sum(bound)) if bound is not None else nf
    out = [{"resnr": atoms[p][0], "resname": atoms[p][1],
            "prot_atom": atoms[p][2], "lig_atom": atoms[l][2],
            "occ": n / nf, "occ_bound": n / max(1, nb)}
           for (p, l), n in hits.items()]
    out.sort(key=lambda r: -r["occ_bound"])
    return out


# --------------------------------------------------------------------------
# 配体聚类
# --------------------------------------------------------------------------
def cluster_ligand(cfg, m: dict, out: Path | None = None, cutoff=0.15):
    """配体构象聚类（GROMOS 法，先按主链叠合再只取配体）。

    gmx cluster 只有一个组选择提示（拟合与 RMSD 用同一组），所以要先用
    trjconv 按主链 rot+trans 叠合、只输出配体，再对这条配体轨迹聚类。
    这里聚的是配体自身构象（cutoff 0.15 nm）；若直接按位姿聚，配体解离后
    会散成上百个碎片簇，没有意义。
    """
    from viz import _gmx
    out = Path(out or m["out"])
    out.mkdir(parents=True, exist_ok=True)
    sfx = m["sfx"]
    lig_gro = out / f"lig_ref{sfx}.gro"
    lig_xtc = out / f"lig_fit{sfx}.xtc"
    sz = out / f"clusters{sfx}.xvg"
    log = out / f"_cluster{sfx}.log"
    if not lig_gro.exists():
        _gmx(cfg, ["trjconv", "-s", m["ref"], "-f", m["ref"], "-n", m["ndx"],
                   "-o", lig_gro, "-dump", "0"],
             cwd=m["job"], stdin_text="LIG\n",
             log=out / f"_lig_ref{sfx}.log")
    if not lig_xtc.exists():
        _gmx(cfg, ["trjconv", "-s", m["ref"], "-f", m["prot"], "-n", m["ndx"],
                   "-o", lig_xtc, "-fit", "rot+trans"],
             cwd=m["job"], stdin_text="Backbone\nLIG\n",
             log=out / f"_lig_fit{sfx}.log")
    if not sz.exists():
        # 不给 -n：让 gmx 用配体结构自己建组（System），否则索引里的组和
        # 24 原子的轨迹对不上；两个提示分别是“聚类组”和“输出组”
        _gmx(cfg, ["cluster", "-s", lig_gro, "-f", lig_xtc,
                   "-g", out / f"cluster{sfx}.log", "-cl",
                   out / f"clusters{sfx}.pdb", "-sz", sz,
                   "-cutoff", cutoff, "-method", "gromos", "-fit", "yes"],
             cwd=m["job"], stdin_text="System\nSystem\n", log=log)
    rows = []
    for ln in sz.read_text().splitlines():
        s = ln.strip()
        if not s or s[0] in "#@":
            continue
        q = s.split()
        if len(q) >= 2:
            rows.append((int(float(q[0])), int(float(q[1]))))
    if not rows:
        return []
    total = sum(n for _, n in rows)
    return [{"cluster": k, "frames": n, "pop": n / total}
            for k, n in rows if n > 0]


# --------------------------------------------------------------------------
# 出图
# --------------------------------------------------------------------------
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = ["Arial", "Microsoft YaHei", "DejaVu Sans"]
    plt.rcParams["axes.linewidth"] = 0.9
    return plt


def plot_contacts(rows, base: Path, title="", top=18, want_pdf=True):
    plt = _plt()
    from make_figures import save
    rows = [r for r in rows if r["occ_bound"] > 0][:top]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(4.6, max(2.8, 0.26 * len(rows) + 1.1)))
    vals = [100 * r["occ_bound"] for r in rows][::-1]
    names = [f"{r['resname'].capitalize()}{r['resnr']}" for r in rows][::-1]
    ax.barh(names, vals, color="#4a7ebb", height=0.66, zorder=2)
    ax.set_xlabel("Bound frames in contact (%)")
    ax.set_xlim(0, 100)
    ax.set_title(title, fontsize=10)
    ax.grid(axis="x", color="#e2e8f0", lw=0.6, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", labelsize=9)
    plt.tight_layout()
    return save(fig, base, want_pdf)


def plot_hbonds(rows, base: Path, title="", top=12, want_pdf=True):
    plt = _plt()
    from make_figures import save
    rows = rows[:top]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(5.2, max(2.4, 0.28 * len(rows) + 1.0)))
    vals = [100 * r["occ_bound"] for r in rows][::-1]
    names = [f"{r['resname'].capitalize()}{r['resnr']} {r['prot_atom']}"
             f"···{r['lig_atom']}" for r in rows][::-1]
    ax.barh(names, vals, color="#2b6cb0", height=0.62, zorder=2)
    ax.set_xlabel("H-bond occupancy in bound frames (%)")
    ax.set_xlim(0, 100)
    ax.set_title(title, fontsize=10)
    ax.grid(axis="x", color="#e2e8f0", lw=0.6, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", labelsize=8.5)
    plt.tight_layout()
    return save(fig, base, want_pdf)


def plot_clusters(rows, base: Path, title="", want_pdf=True):
    plt = _plt()
    from make_figures import save
    rows = rows[:10]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(3.8, 2.8))
    vals = [100 * r["pop"] for r in rows]
    ax.bar([f"C{r['cluster']}" for r in rows], vals, color="#718096",
           width=0.7, zorder=2)
    ax.set_ylabel("Population (%)")
    ax.set_ylim(0, 100)
    ax.set_title(title, fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    plt.tight_layout()
    return save(fig, base, want_pdf)


# --------------------------------------------------------------------------
# 2D 相互作用示意图（RDKit 出骨架，matplotlib 排布残基）
# --------------------------------------------------------------------------
def load_ligand_2d(job: Path, resname="LIG"):
    """读配体 mol2，返回 (mol, 2D 坐标数组, 原子元素列表)。失败返回 None。"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return None
    cand = [job / "01_ligand" / "pose_h.mol2", job / "01_ligand" / "LIG_bcc.mol2"]
    mol = None
    for p in cand:
        if not p.exists():
            continue
        # RDKit 的 C++ 层在 Windows 上打不开中文路径，改成先读成字符串再解析
        block = p.read_text(encoding="utf-8", errors="replace")
        mol = Chem.MolFromMol2Block(block, removeHs=False, sanitize=True)
        if mol is None:
            mol = Chem.MolFromMol2Block(block, removeHs=False, sanitize=False)
        if mol is not None:
            break
    if mol is None:
        return None
    try:
        # 只画重原子骨架：氢会带出一堆悬挂短线，也会把 RDKit 的 2D 排布弄乱
        elems = [a.GetSymbol() for a in mol.GetAtoms()]
        hmap, j = {}, 0
        for i, a in enumerate(mol.GetAtoms()):
            if a.GetAtomicNum() > 1:
                hmap[i] = j
                j += 1
        heavy = Chem.RemoveHs(mol)
        AllChem.Compute2DCoords(heavy)
    except Exception:
        return None
    conf = heavy.GetConformer()
    xy = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                   for i in range(heavy.GetNumAtoms())])
    return heavy, xy, hmap, elems


def map_ligand_names(atoms, resname, lig) -> dict[str, int]:
    """GROMOS 原子名 → mol2 重原子下标；顺序或元素序列对不上就返回 {}。"""
    if lig is None:
        return {}
    _mol, _xy, hmap, elems = lig
    md_lig = [a for a in atoms if a[1] == resname]
    if len(md_lig) != len(elems):
        return {}
    if not all(_element(a[2]) == e.upper() for a, e in zip(md_lig, elems)):
        return {}
    return {a[2]: hmap[k] for k, a in enumerate(md_lig) if k in hmap}


def ligand_anionic_atoms(lig) -> set[int]:
    """找羧酸根：与同一个碳相连、且都是末端（无 H）的两个氧。"""
    if lig is None:
        return set()
    mol = lig[0]
    out = set()
    for c in mol.GetAtoms():
        if c.GetSymbol() != "C":
            continue
        oxy = [n for n in c.GetNeighbors() if n.GetSymbol() == "O"]
        term = [o for o in oxy if o.GetDegree() == 1]
        if len(oxy) >= 2 and len(term) >= 2:
            out.update(o.GetIdx() for o in term)
    return out


def draw_interaction_map(job: Path, m: dict, atoms, hbonds, contacts, bridges,
                         base: Path, title="", want_pdf=True, min_occ=0.02,
                         lig=None):
    """LigPlot 风格径向图；配体缺 mol2 时退化成纯径向图。"""
    plt = _plt()
    from make_figures import save
    lig = lig or load_ligand_2d(job)

    fig, ax = plt.subplots(figsize=(7.0, 5.8))
    if lig is None:
        cx = cy = 0.0
        R = 1.0
        ax.text(0, 0, "LIG", ha="center", va="center", fontsize=14,
                fontweight="bold")
        xy = np.zeros((len(atoms), 2))
    else:
        mol, xy, hmap, elems = lig
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]],
                    color="#1a202c", lw=1.5, solid_capstyle="round", zorder=3)
        for i in range(mol.GetNumAtoms()):
            el = mol.GetAtomWithIdx(i).GetSymbol()
            if el in ("O", "N"):
                ax.text(xy[i, 0], xy[i, 1], el, ha="center", va="center",
                        fontsize=10.5, fontweight="bold", color="#c53030",
                        bbox=dict(boxstyle="circle,pad=0.02", fc="white",
                                  ec="none"), zorder=4)
        cx, cy = xy.mean(0)
        R = 0.5 * max(np.ptp(xy[:, 0]), np.ptp(xy[:, 1])) + 0.6

    # 配体原子名 → 重原子骨架下标；顺序一致且元素序列对得上才挂线
    name2idx = map_ligand_names(atoms, m["res"], lig)

    def merge(rows):
        """同一残基的多个原子对合并，只保留占有率最高的那一对（避免重复计数）。"""
        acc: dict[tuple, dict] = {}
        for r in rows:
            k = (r["resnr"], r["resname"])
            cur = acc.get(k)
            if cur is None or r["occ_bound"] > cur["sum"]:
                acc[k] = {"resnr": k[0], "resname": k[1],
                          "lig_atom": r["lig_atom"],
                          "occ": r["occ_bound"], "sum": r["occ_bound"]}
        return list(acc.values())

    # 每个残基只出现一次：盐桥 > 氢键 > 疏水接触
    chosen: dict[str, tuple] = {}
    for s in merge(bridges):
        lab = f"{s['resname'].capitalize()}{s['resnr']}"
        if s["sum"] >= min_occ:
            chosen[lab] = (lab, min(1.0, s["sum"]), "salt",
                           name2idx.get(s["lig_atom"]))
    for h in merge(hbonds):
        lab = f"{h['resname'].capitalize()}{h['resnr']}"
        if h["sum"] >= min_occ and lab not in chosen:
            chosen[lab] = (lab, min(1.0, h["sum"]), "hbond",
                           name2idx.get(h["lig_atom"]))
    for c in contacts:
        lab = f"{c['resname'].capitalize()}{c['resnr']}"
        if c["occ_bound"] >= 0.15 and lab not in chosen:
            chosen[lab] = (lab, c["occ_bound"], "hydro",
                           name2idx.get(c.get("lig_c_atom", "")))
    items = sorted(chosen.values(), key=lambda t: -t[1])[:12]

    if not items:
        ax.set_xlim(cx - 2 * R, cx + 2 * R)
        ax.set_ylim(cy - 2 * R, cy + 2 * R)
    else:
        # 按挂线原子的方位角排序后均匀铺开，避免标注重叠
        def ang_of(it):
            k = it[3]
            if k is None or k >= len(xy):
                return 0.0
            return float(np.arctan2(xy[k, 1] - cy, xy[k, 0] - cx))
        # 先让标注尽量正对它的目标原子，再按最小间距相互推开（避免连线交叉）
        n = len(items)
        min_sep = 0.92 * 2 * np.pi / max(1, n)
        angs = np.array([ang_of(it) for it in items])
        order = np.argsort(angs)
        items = [items[i] for i in order]
        th = angs[order].astype(float)
        # 从最大空隙之后开始排，避免一圈推下来整体转位（标签跑离它的原子）
        gaps = np.diff(np.r_[th, th[0] + 2 * np.pi])
        start = int(np.argmax(gaps))
        items = items[start:] + items[:start]
        th = np.r_[th[start:], th[:start] + 2 * np.pi]
        for i in range(1, n):
            th[i] = max(th[i], th[i - 1] + min_sep)
        limit = th[0] + 2 * np.pi - min_sep
        if th[-1] > limit and n > 1:
            th -= (th[-1] - limit) * np.linspace(0, 1, n)
        rr_x, rr_y = R + 2.4, R + 1.5
        style = {"hbond": ("#2b6cb0", (0, (4, 2.5)), "-"),
                 "salt": ("#c53030", (0, (5, 2, 1, 2)), "-"),
                 "hydro": ("#718096", (0, (1, 2)), "-")}
        for q, (label, occ, kind, k) in enumerate(items):
            tx = cx + rr_x * np.cos(th[q])
            ty = cy + rr_y * np.sin(th[q])
            color, dash, _ = style[kind]
            if k is not None and k < len(xy):
                ax.plot([tx, xy[k, 0]], [ty, xy[k, 1]], color=color, lw=1.0,
                        ls=dash, zorder=2, alpha=0.95)
            ax.text(tx, ty, f"{label}\n{100 * occ:.0f}%", ha="center",
                    va="center", fontsize=9, color="#1a202c", zorder=5,
                    bbox=dict(boxstyle="round,pad=0.28", fc="white",
                              ec=color, lw=0.9))
        ax.set_xlim(cx - (rr_x + 1.2), cx + (rr_x + 1.2))
        ax.set_ylim(cy - (rr_y + 0.8), cy + (rr_y + 0.8))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=10.5, pad=6)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color="#2b6cb0", ls=(0, (4, 2.5)), lw=1.2,
                      label="H-bond"),
               Line2D([], [], color="#c53030", ls=(0, (5, 2, 1, 2)), lw=1.2,
                      label="Salt bridge"),
               Line2D([], [], color="#718096", ls=(0, (1, 2)), lw=1.4,
                      label="Hydrophobic contact")]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.06),
              ncol=3, frameon=False, fontsize=8.5, handlelength=2.2)
    fig.text(0.5, 0.005, "percentage = fraction of frames with ligand bound",
             ha="center", fontsize=7.5, color="#4a5568")
    plt.tight_layout()
    return save(fig, base, want_pdf)


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def _best_contact_frame(xyz, atoms, masks, cutoff) -> int:
    """接触残基最多（同分时最近）的那一帧，用作手动作图的代表构象。"""
    from scipy.spatial.distance import cdist
    ph, lh = masks["prot_heavy"], masks["lig_heavy"]
    resnr = np.array([atoms[i][0] for i in ph])
    best, bestn, bestd = 0, -1, 1e9
    for f in range(len(xyz)):
        amin = cdist(xyz[f][ph], xyz[f][lh]).min(1)
        n = len(np.unique(resnr[amin <= cutoff]))
        dm = float(amin.min())
        if n > bestn or (n == bestn and dm < bestd):
            best, bestn, bestd = f, n, dm
    return best


def export_pdb(cfg, m: dict, time_ns: float, dst: Path) -> Path:
    """导出指定时刻的 Protein_LIG 复合物 PDB（LigPlot+/PyMOL 都能直接打开）。"""
    from viz import _gmx
    if dst.exists():
        return dst
    _gmx(cfg, ["trjconv", "-s", m["ref"], "-f", m["prot"], "-n", m["ndx"],
               "-o", dst, "-dump", f"{time_ns * 1000:.1f}"],
         cwd=m["job"], stdin_text="Protein_LIG\n",
         log=dst.with_suffix(".log"))
    return dst


def find_ligplot(home: str | None = None) -> Path | None:
    """找 LigPlot+ 的 exe 目录（含 hbplus.exe / ligplot.exe）。

    顺序：显式参数（config.json 的 ligplot_home）→ 环境变量 LIGPLUS_HOME
    → 常见安装位置。任一处不对，把根目录填进 config.json 即可。
    """
    home = home or os.environ.get("LIGPLUS_HOME")
    cands = []
    if home:
        cands += [Path(home), Path(home) / "lib" / "exe_win"]
    cands += [Path(r"D:/LigPlus/LigPlus/lib/exe_win"),
              Path(r"C:/LigPlus/LigPlot+/lib/exe_win"),
              Path(r"C:/Program Files/LigPlus/lib/exe_win")]
    for p in cands:
        if (p / "hbplus.exe").exists() and (p / "ligplot.exe").exists():
            return p
    return None


def run_ligplot(exedir: Path, pdb: Path, outbase: Path, resname="LIG",
                workdir: Path | None = None, dpi=400, want_pdf=True,
                title: str | None = None):
    """用 LigPlot+ 自带的 hbplus.exe + ligplot.exe 出一张经典 LIGPLOT 图。

    LigPlot+ 本体是 Java GUI，但这两个命令行程序可以脱离界面运行：
        hbplus -L -h 2.70 -d 3.90 x.pdb        -> x.hhb
        ligplot x.pdb LIG 136                  -> ligplot.ps
    再用 Ghostscript 转 PNG/PDF 并裁掉白边。

    注意：这样只能画氢键，疏水接触（“睫毛”符号）需要 LigPlot+ 图形界面
    自己算 .nnb，命令行路径拿不到；疏水接触请看本脚本的 interaction_map。
    """
    import shutil
    import subprocess
    work = Path(workdir or (outbase.parent / f"_{outbase.stem}_work"))
    work.mkdir(parents=True, exist_ok=True)
    for exe in ("hbplus.exe", "ligplot.exe"):
        dst = work / exe
        if not dst.exists():
            shutil.copy2(exedir / exe, dst)
    prm = exedir.parent / "params" / "ligplot.prm"
    if prm.exists() and not (work / "ligplot.prm").exists():
        shutil.copy2(prm, work / "ligplot.prm")

    # 配体必须是 HETATM，LigPlot 才认；顺便补上链号，标注才显示 (A)
    lines = []
    for ln in pdb.read_text(encoding="utf-8", errors="replace").splitlines(True):
        if ln.startswith(("ATOM", "HETATM")):
            if ln[17:20] == resname:
                ln = "HETATM" + ln[6:]
        lines.append(ln)
    stem = (title or outbase.stem).replace(" ", "_")
    src = work / f"{stem}.pdb"
    src.write_text("".join(lines), encoding="utf-8")

    # 配体残基号：从 PDB 里取第一个 LIG
    resnum = None
    for ln in lines:
        if ln.startswith("HETATM") and ln[17:20] == resname:
            resnum = int(ln[22:26])
            break
    if resnum is None:
        raise RuntimeError(f"{src.name} 里没有残基 {resname}")

    def _run(args, name):
        p = subprocess.run(args, cwd=work, capture_output=True, text=True,
                           errors="replace")
        (work / f"{name}.log").write_text((p.stdout or "") + "\n" + (p.stderr or ""),
                                          encoding="utf-8")
        return p

    _run([str(work / "hbplus.exe"), "-L", "-h", "2.70", "-d", "3.90", src.name],
         "hbplus")
    _run([str(work / "ligplot.exe"), src.name, resname, str(resnum)], "ligplot")
    ps = work / "ligplot.ps"
    if not ps.exists():
        raise RuntimeError("ligplot.exe 没有生成 ligplot.ps，见 _ligplot_*/ligplot.log")

    gs = shutil.which("gswin64c") or shutil.which("gswin32c") or shutil.which("gs")
    if not gs:
        for c in (r"C:/Program Files/gs/gs10.03.0/bin/gswin64c.exe",
                  r"C:/Program Files/gs/gs10.06.0/bin/gswin64c.exe"):
            if Path(c).exists():
                gs = c
                break
    if not gs:
        raise RuntimeError("没找到 Ghostscript（gswin64c），无法把 .ps 转成图片")
    png = work / "ligplot_raw.png"
    subprocess.run([gs, "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pngalpha",
                    f"-r{dpi}", f"-sOutputFile={png}", str(ps)],
                   capture_output=True)
    if not png.exists():
        raise RuntimeError("Ghostscript 转换失败")
    # 裁白边
    try:
        from PIL import Image, ImageChops
        im = Image.open(png).convert("RGB")
        bg = Image.new("RGB", im.size, (255, 255, 255))
        box = ImageChops.difference(im, bg).getbbox()
        if box:
            pad = 40
            box = (max(0, box[0] - pad), max(0, box[1] - pad),
                   min(im.width, box[2] + pad), min(im.height, box[3] + pad))
            im = im.crop(box)
        im.save(png)
    except Exception:
        pass
    outs = [outbase.with_suffix(".png")]
    outbase.parent.mkdir(parents=True, exist_ok=True)
    png.replace(outs[0]) if png.parent == outbase.parent else shutil.copy2(png, outs[0])
    if want_pdf:
        pdf = outbase.with_suffix(".pdf")
        subprocess.run([gs, "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
                        f"-sOutputFile={pdf}", str(ps)], capture_output=True)
        outs.append(pdf)
    return outs


def run_job(cfg, name: str, rep: int, meta: dict, cutoff=0.45, dt_ns=0.2,
            want_pdf=True, do_map=True, do_ligplot=True) -> dict:
    import viz
    from mdpipe import job_dir
    job = job_dir(name)
    m = viz.ensure_inputs(cfg, job, meta, rep)
    out = job / "08_figures" / "interactions"
    out.mkdir(parents=True, exist_ok=True)
    resname = meta.get("resname", "LIG")
    m["res"] = resname

    atoms = read_gro(m["ref"])
    t, xyz = extract_xyz(cfg, m, dt_ns=dt_ns)
    masks = topology_masks(atoms, resname)

    contacts, bound = contact_table(xyz, atoms, masks, cutoff=cutoff, lig=resname)
    hb = hbond_table(xyz, atoms, masks, resname=resname, bound=bound)
    lig = load_ligand_2d(job) if do_map else None
    neg_names = {n for n, j in map_ligand_names(atoms, resname, lig).items()
                 if j in ligand_anionic_atoms(lig)}
    br = salt_bridges(xyz, atoms, masks, resname=resname, bound=bound,
                      neg_names=neg_names)
    cl = cluster_ligand(cfg, m, out)
    tag = f"_r{rep}" if meta.get("replicates", 1) > 1 else ""

    # 给 LigPlot+ / PyMOL 手动作图留一帧“接触残基最多”的结合构象
    best = _best_contact_frame(xyz, atoms, masks, cutoff)
    export_pdb(cfg, m, t[best], out / f"ligplot{tag}.pdb")
    if do_ligplot:
        exedir = find_ligplot(cfg.get("ligplot_home"))
        if exedir is None:
            print("  没找到 LigPlot+，跳过 LIGPLOT 图（可用 ligplot_home 指定）")
        else:
            try:
                short = name.split("__")[-1].replace("_", " ")
                run_ligplot(exedir, out / f"ligplot{tag}.pdb",
                            out / f"ligplot{tag}", resname=resname,
                            want_pdf=want_pdf,
                            title=f"{short} rep{rep}")
            except Exception as exc:
                print(f"  LigPlot+ 出图失败：{exc}")
    (out / f"contacts{tag}.csv").write_text(
        "resnr,resname,occupancy,occupancy_bound,min_dist_nm,mean_dist_nm\n" +
        "\n".join(f"{r['resnr']},{r['resname']},{r['occ']:.4f},"
                  f"{r['occ_bound']:.4f},{r['min_dist']:.3f},"
                  f"{r['mean_dist']:.3f}" for r in contacts),
        encoding="utf-8")
    (out / f"hbonds{tag}.csv").write_text(
        "resnr,resname,prot_atom,lig_atom,occupancy,occupancy_bound,"
        "mean_dist_nm,mean_angle_deg\n" +
        "\n".join(f"{r['resnr']},{r['resname']},{r['prot_atom']},{r['lig_atom']},"
                  f"{r['occ']:.4f},{r['occ_bound']:.4f},{r['dist']:.3f},"
                  f"{r['angle']:.1f}" for r in hb),
        encoding="utf-8")

    plot_contacts(contacts, out / f"contacts{tag}",
                  title=f"Residue contacts (≤ {cutoff * 10:.1f} Å)", want_pdf=want_pdf)
    plot_hbonds(hb, out / f"hbonds{tag}", title="Protein–ligand hydrogen bonds",
                want_pdf=want_pdf)
    plot_clusters(cl, out / f"clusters{tag}", title="Ligand pose clusters",
                  want_pdf=want_pdf)
    if do_map:
        draw_interaction_map(job, m, atoms, hb, contacts, br,
                             out / f"interaction_map{tag}", lig=lig,
                             title=f"{name.split('__')[-1].replace('_', ' ')}"
                                   f"  ·  {t[0]:.0f}–{t[-1]:.0f} ns",
                             want_pdf=want_pdf)

    summary = {
        "job": name, "rep": rep, "resname": resname,
        "window_ns": [float(t[0]), float(t[-1])], "frames": int(len(t)),
        "bound_fraction": float(np.mean(bound)),
        "cutoff_nm": cutoff,
        "top_contacts": contacts[:12],
        "hbonds": [h for h in hb if h["occ_bound"] >= 0.05],
        "salt_bridges": [s for s in br if s["occ_bound"] >= 0.05],
        "clusters": cl[:6],
    }
    (out / f"interactions{tag}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description="配体-蛋白相互作用分析")
    ap.add_argument("--job", nargs="+", required=True, help="流水线作业名")
    ap.add_argument("--replicates", nargs="+", type=int)
    ap.add_argument("--cutoff", type=float, default=0.45, help="接触截断 nm")
    ap.add_argument("--dt", type=float, default=0.2, help="采样间隔 ns")
    ap.add_argument("--no-map", action="store_true", help="不画 2D 示意图")
    ap.add_argument("--no-ligplot", action="store_true",
                    help="不调用 LigPlot+（默认检测到就调用）")
    ap.add_argument("--no-pdf", action="store_true")
    a = ap.parse_args(argv)

    import mdpipe
    cfg = mdpipe.config()
    for name in a.job:
        meta = mdpipe.job_meta(name)
        nrep = int(meta.get("replicates", 1))
        reps = a.replicates or list(range(1, nrep + 1))
        for rep in reps:
            print(f"[{name}] 重复 {rep} …", flush=True)
            s = run_job(cfg, name, rep, meta, cutoff=a.cutoff, dt_ns=a.dt,
                        want_pdf=not a.no_pdf, do_map=not a.no_map,
                        do_ligplot=not a.no_ligplot)
            hb = s["hbonds"][0] if s["hbonds"] else None
            ncon = sum(1 for c in s["top_contacts"] if c["occ_bound"] >= 0.25)
            print(f"  结合帧占比 {100 * s['bound_fraction']:.0f}%；"
                  f"接触残基 {ncon} 个（结合帧内 ≥25%）；"
                  f"氢键 {len(s['hbonds'])} 对"
                  + (f"，最强 {hb['resname']}{hb['resnr']} {hb['prot_atom']}"
                     f"···{hb['lig_atom']} {100 * hb['occ_bound']:.0f}%" if hb else ""),
                  flush=True)


if __name__ == "__main__":
    main()
