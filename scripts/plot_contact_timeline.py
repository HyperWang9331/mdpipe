#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""各体系的"配体在口袋内"时间线：接触原子数 vs 时间，两个重复画在同一小图里。

接触原子数 = 配体 0.35 nm 内的蛋白原子数（gmx mindist -on），取 5 ns 滑动平均；
500 是"在口袋里"的判定线（FABP5 口袋满员约 900–1350 个原子）。

    python scripts/plot_contact_timeline.py
"""

from __future__ import annotations

import glob
import os
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
JOBS = ROOT / "pipeline" / "jobs"
CUT = 500

# 画图顺序与标签；名字用于匹配作业目录
ORDER = [("D_Galacturonic_Acid", "D-GalA", -6.2), ("D_Galactose", "D-Gal", -5.5),
         ("D_Glucuronic_Acid", "D-GlcA", -5.5), ("D_Mannose", "D-Man", -5.1),
         ("D_Glucose", "D-Glc", -5.1), ("L_Arabinose", "L-Ara", -4.6),
         ("L_Rhamnose", "L-Rha", -5.1)]
COLOR = {"D-GalA": "#2b6cb0", "D-Gal": "#2b6cb0", "D-GlcA": "#2b6cb0",
         "D-Man": "#c05621", "D-Glc": "#c05621", "L-Ara": "#c05621",
         "L-Rha": "#c05621"}


def load(job: str, rep: int):
    d = glob.glob(str(JOBS / f"*{job}*"))
    if not d:
        return None, None
    p = pathlib.Path(d[0]) / f"07_analysis_r{rep}" / "contacts.xvg"
    if not p.exists():
        return None, None
    x = np.loadtxt(p, comments=["#", "@"])
    t, c = x[:, 0] / 1000.0, x[:, 1]
    k = max(1, int(5.0 / (t[1] - t[0])))
    return t[k - 1:], np.convolve(c, np.ones(k) / k, mode="valid")


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    fig, axes = plt.subplots(2, 4, figsize=(13.0, 6.0), dpi=300,
                             sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, (job, tag, dock) in zip(axes, ORDER):
        got = False
        for r in (1, 2):
            t, sm = load(job, r)
            if t is None:
                continue
            got = True
            ax.plot(t, sm, lw=1.5 if r == 1 else 1.0, color=COLOR[tag],
                    ls="-" if r == 1 else "--", alpha=0.95 if r == 1 else 0.65,
                    label=f"replicate {r}")
        ax.axhspan(0, CUT, color="#f2f4f7", zorder=0)
        ax.axhline(CUT, color="#a0aec0", lw=0.8, ls=":", zorder=1)
        ax.set_title(f"{tag}   (docking {dock:+.1f})", fontsize=9.5)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 1450)
        ax.set_xticks([0, 50, 100])
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        if not got:
            ax.text(0.5, 0.55, "MD in progress", transform=ax.transAxes,
                    ha="center", va="center", fontsize=9, color="#a0aec0")
        elif ax is axes[0]:
            ax.legend(frameon=False, fontsize=7.5, loc="lower left")

    for ax in axes[len(ORDER):]:
        ax.axis("off")
    for ax in axes[4:7]:
        ax.set_xlabel("Time (ns)", fontsize=9)
    for ax in (axes[0], axes[4]):
        ax.set_ylabel("Protein atoms within\n0.35 nm of ligand", fontsize=9)

    fig.text(0.5, 0.012,
             "Shaded band = ligand out of the canonical pocket (<500 contact atoms); "
             "curves are 5-ns running means",
             ha="center", fontsize=8, color="#718096")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out = ROOT / "pipeline" / "figures" / "contact_timeline.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    print(f"写出 {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
