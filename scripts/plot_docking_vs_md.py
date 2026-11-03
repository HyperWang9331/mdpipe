#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""对接打分 vs MD 口袋停留率。

左边喂 Vina 的 output.pdbqt，右边读流水线分析出来的 contacts.xvg，
按配体名对上，画一张「静态打分能不能预测动态稳定性」的散点图。

    python scripts/plot_docking_vs_md.py \
        --docking "C:/.../对接目录" --out pipeline/figures/docking_vs_md.png
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = pathlib.Path(__file__).resolve().parent.parent
JOBS = ROOT / "pipeline" / "jobs"

# 作业名里的标签 ↔ 对接目录名
ALIAS = {
    "D_Glucose": "D-Glc", "D_Galactose": "D-Gal", "D_Mannose": "D-Man",
    "L_Arabinose": "L-Ara", "L_Rhamnose": "L-Rha",
    "D_Glucuronic_Acid": "D-GlcA", "D_Galacturonic_Acid": "D-GalA",
}
POCKET_CUT = 500          # 0.35 nm 内蛋白原子数 ≥ 500 视为"在口袋里"


def docking_scores(root: pathlib.Path) -> dict[str, float]:
    out = {}
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        f = d / "output.pdbqt"
        if not f.exists():
            continue
        m = re.search(r"REMARK VINA RESULT:\s*(-?\d+\.\d+)", f.read_text(errors="replace"))
        if not m:
            continue
        key = d.name.split("+")[0]
        out[ALIAS.get(key, key)] = float(m.group(1))
    return out


def md_occupancy() -> dict[str, list[tuple[int, float, float]]]:
    """{标签: [(重复号, 在口袋帧占比, 后50ns平均接触原子数), ...]}"""
    res: dict[str, list] = {}
    for job in JOBS.iterdir():
        if not job.is_dir():
            continue
        for key in sorted(ALIAS, key=len, reverse=True):
            if key not in job.name:
                continue
            tag = ALIAS[key]
            reps = []
            for r in (1, 2):
                p = job / f"07_analysis_r{r}" / "contacts.xvg"
                if not p.exists():
                    continue
                d = np.loadtxt(p, comments=["#", "@"])
                t, n = d[:, 0], d[:, 1]
                late = n[t >= 50000]
                if late.size == 0:
                    continue
                reps.append((r, 100.0 * float((late >= POCKET_CUT).mean()), float(late.mean())))
            if reps:
                res[tag] = reps
                break
    return res


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description="对接打分 vs MD 口袋停留率")
    ap.add_argument("--docking", required=True, help="对接结果根目录（每个子目录一个 output.pdbqt）")
    ap.add_argument("--out", default=str(ROOT / "pipeline" / "figures" / "docking_vs_md.png"))
    a = ap.parse_args(argv)

    dock = docking_scores(pathlib.Path(a.docking))
    md = md_occupancy()
    done = {k: v for k, v in md.items() if k in dock}
    pending = sorted(k for k in dock if k not in md)

    fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=300)
    ax.axhspan(-8, 50, color="#f7fafc", zorder=0)
    ax.axhline(50, color="#cbd5e0", lw=0.9, ls="--", zorder=0)
    ax.text(0.015, 0.055, "left the pocket in the late 50 ns", transform=ax.transAxes,
            fontsize=8, color="#718096")

    # 打平分的体系在横轴上错开一点；标签位置逐体系手调，避免叠在一起
    JITTER = {"D-Gal": -0.035, "D-GlcA": 0.035}
    # (横轴偏移, 纵轴偏移) 均为数据单位，后两项为对齐方式
    LABEL = {"D-GalA": (-0.02, -5, "center", "top"),
             "D-Gal": (-0.06, 0, "left", "center"),
             "D-GlcA": (0.06, 0, "right", "center"),
             "D-Glc": (0.025, 3, "left", "center")}

    for tag, reps in done.items():
        x = dock[tag] + JITTER.get(tag, 0.0)
        ys = [r[1] for r in reps]
        color = "#2b6cb0" if tag.endswith("A") else "#c05621"
        if len(ys) == 2:
            ax.plot([x, x], ys, color=color, lw=1.0, alpha=0.5, zorder=2)
        # 先画空心的重复 2、再画实心的重复 1：两个值相同时看见的是一个套环，表示"两次一致"
        for rep, occ, _ in sorted(reps, key=lambda q: q[0], reverse=True):
            filled = rep == 1
            ax.scatter([x], [occ], s=70 if filled else 150,
                       facecolor=color if filled else "white",
                       edgecolor=color, linewidth=1.2, zorder=3)
        dx, dy, ha, va = LABEL.get(tag, (0.0, 8, "center", "bottom"))
        ax.annotate(tag, (x, max(ys)), xytext=(x + dx, max(ys) + dy), textcoords="data",
                    fontsize=8.5, color="#2d3748", ha=ha, va=va)
        if len(ys) == 2 and abs(ys[0] - ys[1]) >= 15:
            # 数值标在空心（重复 2）那一点旁边，避免被误读成实心点的值
            ax.annotate(f"{ys[1]:.0f}%", (x, ys[1]), xytext=(12, -3),
                        textcoords="offset points", fontsize=8, color=color, va="center")

    if pending:
        ax.text(0.985, 0.035, "MD still running: " + ", ".join(pending),
                transform=ax.transAxes, ha="right", fontsize=7.5, color="#a0aec0")

    ax.set_xlabel(r"Docking $\Delta G$ (kcal mol$^{-1}$)", fontsize=10)
    ax.set_ylabel("Frames in canonical pocket (%)", fontsize=10)
    ax.set_ylim(-8, 108)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.invert_xaxis()
    ax.grid(axis="y", color="#e2e8f0", lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    handles = [Line2D([], [], marker="o", ls="", ms=6, color="#4a5568", label="replicate 1"),
               Line2D([], [], marker="o", ls="", ms=9, mfc="white", mec="#4a5568",
                      label="replicate 2 (ring = both agree)")]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="center right",
              bbox_to_anchor=(1.0, 0.72), handletextpad=0.5, labelspacing=0.7)

    fig.tight_layout()
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, bbox_inches="tight", facecolor="white")
    print(f"写出 {a.out}")
    print("\n体系      对接ΔG   在口袋帧占比（重复1 / 重复2）")
    for tag in sorted(done, key=lambda k: dock[k]):
        reps = done[tag]
        print(f"  {tag:7s} {dock[tag]:6.1f}   " + " / ".join(f"{r[1]:.0f}%" for r in reps))
    if pending:
        print("  未完成：" + ", ".join(f"{t} ({dock[t]:.1f})" for t in pending))
    return 0


if __name__ == "__main__":
    sys.exit(main())
