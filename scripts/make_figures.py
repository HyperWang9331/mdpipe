#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_figures.py —— 把 MD 分析结果画成论文用的图

单个体系：六联图（骨架 RMSD / 配体 RMSD / Rg / SASA / 氢键 / RMSF），多次重复叠画。
多个体系：同一套面板画在一张图上，每个体系一条均值曲线 + 重复范围带，用来横向比较。

用法
----
    python scripts/make_figures.py --job Gal_1          # pipeline/jobs/Gal_1
    python scripts/make_figures.py --dir path/to/job    # 任意目录
    python scripts/make_figures.py --all                # 队列里所有作业
    python scripts/make_figures.py --compare Gal_1 Glc_1 Man_1
    python scripts/make_figures.py --job Gal_1 --no-pdf # 只要 PNG

输入是 mdpipe 分析步骤产出的 xvg：rmsd_backbone / rmsd_ligand / rg / sasa / hbond / rmsf。
输出：<作业目录>/08_figures/<作业名>_MD_analysis.{png,pdf,svg} 与统计 CSV。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib import rcParams                      # noqa: E402

rcParams["font.family"] = ["Arial", "Microsoft YaHei", "DejaVu Sans"]
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 9.5
rcParams["axes.linewidth"] = 0.9
rcParams["xtick.direction"] = "out"
rcParams["ytick.direction"] = "out"
rcParams["axes.grid"] = False
SERIF = "Times New Roman"

RED, BLUE, GREY = "#FF0000", "#0000FF", "#808080"
PALETTE = ["#FF0000", "#0000FF", "#808080", "#008000", "#FF7F00",
           "#800080", "#8B4513", "#00CED1", "#FF1493", "#2F4F4F"]

# 面板定义：(标题, 文件名模板, y 轴, 取第几列)
PANELS = [
    ("Backbone RMSD", "rmsd_backbone_{r}.xvg", "Backbone RMSD (nm)", 1),
    ("Ligand RMSD", "rmsd_ligand_{r}.xvg", "Ligand RMSD (nm)", 1),
    ("Rg", "rg_{r}.xvg", "Radius of gyration (nm)", 1),
    ("SASA", "sasa_{r}.xvg", "SASA (nm$^2$)", 1),
    # gmx hbond-legacy -num 的第 2 列 "Hydrogen bonds" 才是氢键条数，第 3 列是
    # "Pairs within 0.35 nm"（距离判据内的原子对数，恒 ≥ 氢键数），别取错。
    ("H-bonds", "hbond_{r}.xvg", "Number of H-bonds", 1),
    ("RMSF", "rmsf_{r}.xvg", "RMSF (nm)", 1),
]


def load_xvg(path: Path, col: int = 1):
    """读 xvg，col 从 0 开始（0 是时间列）。col=-1 表示取最后一列。"""
    if not path or not Path(path).exists():
        return None, None
    xs, ys = [], []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s or s[0] in "#@":
            continue
        q = s.split()
        c = len(q) - 1 if col == -1 else col
        if len(q) <= c:
            continue
        try:
            xs.append(float(q[0]))
            ys.append(float(q[c]))
        except ValueError:
            continue
    if not xs:
        return None, None
    return np.array(xs), np.array(ys)


def analysis_dirs(job: Path) -> list[Path]:
    """找出 07_analysis、07_analysis_r1、07_analysis_r2 ... 按重复编号排序。"""
    dirs = [d for d in job.glob("07_analysis*") if d.is_dir()]
    return sorted(dirs, key=lambda p: (p.name != "07_analysis", p.name))


def find_one(d: Path, tpl: str):
    """先找带重复编号的（rmsd_backbone_1.xvg），再找不带编号的旧文件名。"""
    pats = [tpl, tpl.replace("rmsd_ligand", "rmsd_lig")]   # 旧版流水线叫 rmsd_lig
    for t in pats:
        for pat in (t.format(r="*"), t.replace("_{r}", "")):
            cand = sorted(d.glob(pat))
            if cand:
                return cand[0]
    return None


def collect(job: Path):
    """返回 {文件名模板: [(x, y), ...]}，一个元素代表一次重复。"""
    data = {tpl: [] for _, tpl, _, _ in PANELS}
    for d in analysis_dirs(job):
        for _, tpl, _, col in PANELS:
            p = find_one(d, tpl)
            if p is None:
                continue
            x, y = load_xvg(p, col)
            if y is not None:
                data[tpl].append((x, y))
    return data


def style_panel(ax, letter: str):
    ax.text(-0.18, 1.06, letter, transform=ax.transAxes, fontsize=14,
            fontweight="bold", fontfamily=SERIF, va="bottom", ha="left")


def save(fig, base: Path, want_pdf=True):
    outs = [base.with_suffix(".png")]
    if want_pdf:
        outs += [base.with_suffix(".pdf"), base.with_suffix(".svg")]
    for p in outs:
        fig.savefig(p, dpi=500 if p.suffix == ".png" else None,
                    bbox_inches="tight", facecolor="white", pad_inches=0.15)
    return outs


def plot_job(job: Path, name: str, want_pdf=True):
    data = collect(job)
    if not any(data.values()):
        print(f"{name}: 没有找到 xvg，先跑分析（mdpipe 的分析步骤）")
        return None

    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.8))
    stats = []
    flat = [axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 0], axes[1, 1], axes[1, 2]]

    for i, ((title, tpl, ylab, _), ax) in enumerate(zip(PANELS, flat)):
        curves = data[tpl]
        if not curves:
            ax.set_axis_off()
            continue
        multi = len(curves) > 1
        for k, (x, y) in enumerate(curves):
            color = [RED, BLUE, GREY, "#008000"][k % 4]
            ax.plot(x, y, color=color, lw=0.6,
                    label=f"Replicate {k + 1}" if multi else None)
            tail = y[int(len(y) * 0.5):]
            stats.append((title, k + 1, y.mean(), y.std(), tail.mean(), tail.std()))
        ax.set_xlabel("Residue index" if title == "RMSF" else "Time (ns)")
        ax.set_ylabel(ylab)
        style_panel(ax, "ABCDEF"[i])
        if multi:
            ax.legend(frameon=False, fontsize=7.5, loc="upper right",
                      handlelength=1.6, borderaxespad=0.3)

    plt.subplots_adjust(wspace=0.34, hspace=0.42)
    outdir = job / "08_figures"
    outdir.mkdir(parents=True, exist_ok=True)
    outs = save(fig, outdir / f"{name}_MD_analysis", want_pdf)
    plt.close(fig)

    csv = outdir / f"{name}_MD统计.csv"
    with open(csv, "w", encoding="utf-8-sig") as fh:
        fh.write("Metric,Replicate,Mean,SD,Tail50% mean,Tail50% SD\n")
        for r in stats:
            fh.write(f"{r[0]},{r[1]},{r[2]:.4f},{r[3]:.4f},{r[4]:.4f},{r[5]:.4f}\n")
    print(f"{name}: {', '.join(p.name for p in outs)}  + {csv.name}")
    for tag in ("Backbone RMSD", "Ligand RMSD"):
        vals = [s for s in stats if s[0] == tag]
        if vals:
            print(f"    {tag}: " + ", ".join(
                f"重复{int(v[1])} {v[2]:.3f}±{v[3]:.3f} nm" for v in vals))
    return outs[0]


def plot_compare(jobs: list[Path], names: list[str], outdir: Path, want_pdf=True):
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.8))
    flat = [axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 0], axes[1, 1], axes[1, 2]]
    plotted = False

    for i, ((title, tpl, ylab, _), ax) in enumerate(zip(PANELS, flat)):
        for k, (job, nm) in enumerate(zip(jobs, names)):
            data = collect(job).get(tpl, [])
            if not data:
                continue
            n = min(len(y) for _, y in data)
            x = data[0][0][:n]
            M = np.vstack([y[:n] for _, y in data])
            color = PALETTE[k % len(PALETTE)]
            if len(data) > 1:
                ax.fill_between(x, M.min(0), M.max(0), color=color, alpha=0.15, lw=0)
            ax.plot(x, M.mean(0), color=color, lw=1.0, label=nm)
            plotted = True
        ax.set_xlabel("Residue index" if title == "RMSF" else "Time (ns)")
        ax.set_ylabel(ylab)
        style_panel(ax, "ABCDEF"[i])
        if ax.get_legend_handles_labels()[0]:
            ax.legend(frameon=False, fontsize=7, loc="upper right", handlelength=1.6)

    if not plotted:
        print("没有找到可画的数据")
        plt.close(fig)
        return None
    plt.subplots_adjust(wspace=0.34, hspace=0.42)
    outdir.mkdir(parents=True, exist_ok=True)
    outs = save(fig, outdir / "MD_compare", want_pdf)
    plt.close(fig)
    print("对比图: " + ", ".join(p.name for p in outs))
    return outs[0]


def main():
    ap = argparse.ArgumentParser(description="MD 分析出图")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--job", help="队列里的作业名")
    g.add_argument("--dir", help="任意作业目录")
    g.add_argument("--all", action="store_true", help="队列里所有作业")
    g.add_argument("--compare", nargs="+", metavar="作业", help="多体系对比图")
    ap.add_argument("--name", help="图里显示的体系名，默认用作业名")
    ap.add_argument("--outdir", default=None, help="对比图输出目录")
    ap.add_argument("--no-pdf", action="store_true", help="只出 PNG")
    a = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import mdpipe

    if a.compare:
        jobs = [mdpipe.job_dir(n) for n in a.compare]
        missing = [n for n, j in zip(a.compare, jobs) if not j.exists()]
        if missing:
            sys.exit(f"找不到作业目录：{missing}")
        outdir = Path(a.outdir) if a.outdir else mdpipe.PIPE / "figures"
        plot_compare(jobs, a.compare, outdir, not a.no_pdf)
        return

    if a.dir:
        job = Path(a.dir)
        name = a.name or job.name
        if not job.exists():
            sys.exit(f"目录不存在：{job}")
        plot_job(job, name, not a.no_pdf)
        return

    names = mdpipe.queue()["jobs"] if a.all else [a.job]
    for n in names:
        job = mdpipe.job_dir(n)
        if not job.exists():
            print(f"跳过 {n}：目录不存在")
            continue
        plot_job(job, a.name or n, not a.no_pdf)


if __name__ == "__main__":
    main()
