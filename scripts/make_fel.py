#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_fel.py —— PCA 自由能景观（FEL）出图

数据来自 mdpipe 的 pca 步骤（每个重复一个 pca/ 目录）：
    eigenval.xvg   特征值谱
    proj.xvg       PC1、PC2 随时间的投影
    2dproj.xvg     PC1–PC2 二维投影（算 FEL 用）

出两张图：
    <作业名>_PCA_FEL.png/pdf/svg   四联图：特征值谱、投影随时间、每重复的 2D 景观
    <作业名>_FEL_3D.png/pdf/svg    3D 曲面 + 地面二维映射（论文常用的那张）

多体系对比：
    <outdir>/FEL_compare_3D.png/pdf/svg   1×N 3D 曲面，共用一个色标

约定与原来的出图代码一致：Spectral 色板、45×45 网格（帧少时自动减格）、σ=0.9 平滑、
ΔG = −kBT ln(P/Pmax)、kT = 0.0083144621 × 300 kJ/mol、上限 15 kJ/mol。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                          # noqa: E402
from matplotlib import rcParams                          # noqa: E402
from matplotlib.colors import LinearSegmentedColormap    # noqa: E402

try:                                                     # 平滑可选，缺 scipy 就跳过
    from scipy.ndimage import gaussian_filter
except Exception:                                        # noqa: BLE001
    gaussian_filter = None

rcParams["font.family"] = ["Arial", "Microsoft YaHei", "DejaVu Sans"]
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 9.5
rcParams["axes.linewidth"] = 0.9
rcParams["xtick.direction"] = "out"
rcParams["ytick.direction"] = "out"
rcParams["axes.grid"] = False
SERIF = "Times New Roman"

RED, BLUE = "#FF0000", "#0000FF"
KB, T = 0.0083144621, 300.0
NBINS, VMAX = 45, 15.0
SPECTRAL = LinearSegmentedColormap.from_list(
    "spectral_ref",
    ["#5e4fa2", "#3288bd", "#66c2a5", "#abdda4", "#e6f598",
     "#fee08b", "#fdae61", "#f46d43", "#d53e4f"], N=256)


def load_xvg(path):
    rows = []
    if not path.exists():
        return None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s or s[0] in "#@&":
            continue
        q = s.split()
        try:
            rows.append([float(x) for x in q])
        except ValueError:
            continue
    return np.array(rows) if rows else None


def load_proj(path):
    """anaeig -proj 每条特征向量一块，用 '&' 分隔；返回 [PC1 数组, PC2 数组]。"""
    blocks, cur = [], []
    if not path.exists():
        return None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if s.startswith("&"):
            if cur:
                blocks.append(np.array(cur))
                cur = []
            continue
        if not s or s[0] in "#@":
            continue
        q = s.split()
        try:
            cur.append([float(q[0]), float(q[1])])
        except (ValueError, IndexError):
            continue
    if cur:
        blocks.append(np.array(cur))
    return blocks or None


def analysis_dirs(job: Path):
    dirs = [d for d in job.glob("07_analysis*") if d.is_dir() and (d / "pca").is_dir()]
    return sorted(dirs, key=lambda p: (p.name != "07_analysis", p.name))


def collect(job: Path):
    """每个重复一份：{eig, proj, proj2d}。"""
    reps = []
    for d in analysis_dirs(job):
        p = d / "pca"
        reps.append({
            "dir": d.name,
            "eig": load_xvg(p / "eigenval.xvg"),
            "proj": load_proj(p / "proj.xvg"),
            "proj2d": load_xvg(p / "2dproj.xvg"),
        })
    return reps


def fel_grid(d, nbins=None):
    """二维投影 → 自由能面。空网格按上限填满，和原脚本一致。

    长轨迹用 45×45；帧数少于 4×NBINS 时自动减网格，否则满屏都是空白格被顶到
    ΔG 上限，曲面会变成一片"天花板"，看不出盆地。
    """
    if nbins is None:
        nbins = NBINS if len(d) >= 4 * NBINS else max(12, int(np.sqrt(len(d)) * 2))
    H, xe, ye = np.histogram2d(d[:, 0], d[:, 1], bins=nbins)
    if gaussian_filter is not None:
        H = gaussian_filter(H, sigma=0.9)
    P = H / H.sum() if H.sum() else H
    with np.errstate(divide="ignore", invalid="ignore"):
        G = -KB * T * np.log(P / P.max())
    G[~np.isfinite(G)] = VMAX
    G = np.clip(G, 0, VMAX)
    xc = 0.5 * (xe[:-1] + xe[1:])
    yc = 0.5 * (ye[:-1] + ye[1:])
    X, Y = np.meshgrid(xc, yc)
    return X, Y, G


def tag(ax, s, x=-0.175, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=14, fontweight="bold",
            fontfamily=SERIF, va="bottom", ha="left")


def save(fig, base: Path, want_pdf=True):
    outs = [base.with_suffix(".png")]
    if want_pdf:
        outs += [base.with_suffix(".pdf"), base.with_suffix(".svg")]
    for p in outs:
        fig.savefig(p, dpi=500 if p.suffix == ".png" else None,
                    bbox_inches="tight", facecolor="white", pad_inches=0.15)
    return outs


def plot_job(job: Path, name: str, want_pdf=True):
    reps = [r for r in collect(job) if r["proj2d"] is not None]
    if not reps:
        print(f"{name}: 没有 PCA 数据（轨迹太短或没跑 pca 步骤）")
        return None

    n = len(reps)
    if n == 1:
        fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.3))
        axes = np.atleast_2d(axes)
        cells = [axes[0, 0], axes[0, 1], axes[0, 2]]
        letters = "ABC"
    else:
        fig, axes = plt.subplots(2, 2, figsize=(9.8, 8.2))
        cells = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]]
        letters = "ABCD"

    colors = [RED, BLUE, "#808080", "#008000"]

    # A 特征值谱
    ax = cells[0]
    for k, r in enumerate(reps):
        if r["eig"] is None:
            continue
        m = min(len(r["eig"]), 30)
        ax.plot(range(1, m + 1), r["eig"][:m, 1], color=colors[k % 4], lw=1.0,
                marker="os"[k % 2], ms=3.0, markerfacecolor="white",
                markeredgewidth=1.0, label=f"Replicate {k + 1}")
    ax.set_xlabel("Eigenvector index")
    ax.set_ylabel("Eigenvalue (nm$^2$)")
    ax.legend(frameon=False, fontsize=8)
    tag(ax, letters[0])

    # B 投影随时间
    ax = cells[1]
    for k, r in enumerate(reps):
        if not r["proj"]:
            continue
        c = colors[k % 4]
        b = r["proj"]
        ax.plot(b[0][:, 0] / 1000.0, b[0][:, 1], color=c, lw=0.65,
                label=f"Replicate {k + 1}  PC1")
        if len(b) >= 2:
            ax.plot(b[1][:, 0] / 1000.0, b[1][:, 1], color=c, lw=0.65, ls="--",
                    alpha=0.85, label=f"Replicate {k + 1}  PC2")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Projection (nm)")
    ax.legend(frameon=False, fontsize=7.5)
    tag(ax, letters[1])

    # C 起：每个重复一张 2D 景观
    for k, r in enumerate(reps):
        if 2 + k >= len(cells):
            break
        ax = cells[2 + k]
        X, Y, G = fel_grid(r["proj2d"])
        im = ax.pcolormesh(X, Y, G, cmap=SPECTRAL, shading="auto",
                           vmin=0, vmax=VMAX)
        cb = plt.colorbar(im, ax=ax, pad=0.02)
        cb.set_label("Free energy (kJ/mol)", fontsize=8.5)
        cb.ax.tick_params(labelsize=8)
        ax.set_xlabel("PC1 (nm)")
        ax.set_ylabel("PC2 (nm)")
        ax.text(0.97, 0.96, f"Replicate {k + 1}", transform=ax.transAxes,
                fontsize=8.5, ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.75",
                          lw=0.6, alpha=0.9))
        tag(ax, letters[2 + k])

    plt.subplots_adjust(wspace=0.34, hspace=0.40)
    outdir = job / "08_figures"
    outdir.mkdir(parents=True, exist_ok=True)
    outs = save(fig, outdir / f"{name}_PCA_FEL", want_pdf)
    plt.close(fig)

    # 3D 曲面（所有重复的点合成一张景观）
    allpts = np.vstack([r["proj2d"] for r in reps])
    fig = plt.figure(figsize=(6.6, 5.4))
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    X, Y, G = fel_grid(allpts)
    surf = ax.plot_surface(X, Y, G, cmap=SPECTRAL, vmin=0, vmax=VMAX,
                           linewidth=0, antialiased=True,
                           rcount=G.shape[0], ccount=G.shape[1], shade=False)
    ax.contourf(X, Y, G, zdir="z", offset=0.0, levels=18, cmap=SPECTRAL,
                vmin=0, vmax=VMAX, antialiased=True)
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.set_zlim(0, VMAX)
    ax.set_box_aspect((1, 1, 0.70), zoom=1.02)
    ax.set_xlabel("PC1 (nm)", fontsize=9, labelpad=-4)
    ax.set_ylabel("PC2 (nm)", fontsize=9, labelpad=-4)
    ax.view_init(elev=26, azim=-58)
    ax.tick_params(labelsize=7.5, pad=-1.5)
    cb = fig.colorbar(surf, shrink=0.62, pad=0.08)
    cb.set_label("Free energy (kJ/mol)", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    outs3 = save(fig, outdir / f"{name}_FEL_3D", want_pdf)
    plt.close(fig)

    print(f"{name}: {', '.join(p.name for p in outs + outs3)}")
    for k, r in enumerate(reps):
        d = r["proj2d"]
        print(f"    重复{k + 1}: {len(d)} 帧 | PC1 {d[:, 0].min():.2f}~{d[:, 0].max():.2f} nm"
              f" | PC2 {d[:, 1].min():.2f}~{d[:, 1].max():.2f} nm")
    return outs[0]


def plot_compare(jobs: list[Path], names: list[str], outdir: Path, want_pdf=True):
    """1×N 三维曲面，共用一个色标（对照图18的做法）。"""
    data = []
    for job, nm in zip(jobs, names):
        reps = [r for r in collect(job) if r["proj2d"] is not None]
        if not reps:
            print(f"跳过 {nm}：没有 PCA 数据")
            continue
        allpts = np.vstack([r["proj2d"] for r in reps])
        data.append((nm, fel_grid(allpts)))
    if not data:
        return None

    n = len(data)
    fig = plt.figure(figsize=(4.0 * n + 1.2, 5.0))
    last = None
    for k, (nm, (X, Y, G)) in enumerate(data):
        ax = fig.add_subplot(1, n, k + 1, projection="3d")
        last = ax.plot_surface(X, Y, G, cmap=SPECTRAL, vmin=0, vmax=VMAX,
                               linewidth=0, antialiased=True,
                               rcount=G.shape[0], ccount=G.shape[1], shade=False)
        ax.contourf(X, Y, G, zdir="z", offset=0.0, levels=18, cmap=SPECTRAL,
                    vmin=0, vmax=VMAX, antialiased=True)
        ax.set_xlim(X.min(), X.max())
        ax.set_ylim(Y.min(), Y.max())
        ax.set_zlim(0, VMAX)
        ax.set_box_aspect((1, 1, 0.70), zoom=1.02)
        ax.set_xlabel("PC1 (nm)", fontsize=9, labelpad=-4)
        ax.set_ylabel("PC2 (nm)", fontsize=9, labelpad=-4)
        ax.view_init(elev=26, azim=-58)
        ax.tick_params(labelsize=7.5, pad=-1.5)
        ax.text2D(-0.09, 1.00, "ABCDEFGH"[k], transform=ax.transAxes, fontsize=14,
                  fontweight="bold", fontfamily=SERIF, va="bottom", ha="left")
        ax.set_title(nm, fontsize=9.5, pad=0)

    fig.subplots_adjust(left=0.005, right=0.845, bottom=0.01, top=0.99, wspace=0.10)
    cax = fig.add_axes([0.895, 0.26, 0.021, 0.50])
    cb = fig.colorbar(last, cax=cax)
    cb.set_label("Free energy (kJ/mol)", fontsize=10, labelpad=3)
    cb.ax.tick_params(labelsize=8.5, pad=2)

    outdir.mkdir(parents=True, exist_ok=True)
    outs = save(fig, outdir / "FEL_compare_3D", want_pdf)
    plt.close(fig)
    print("FEL 对比图: " + ", ".join(p.name for p in outs))
    return outs[0]


def main():
    ap = argparse.ArgumentParser(description="PCA 自由能景观出图")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--job", help="队列里的作业名")
    g.add_argument("--dir", help="任意作业目录")
    g.add_argument("--compare", nargs="+", metavar="作业", help="多体系 FEL 对比")
    ap.add_argument("--name", help="图里显示的名字")
    ap.add_argument("--outdir", help="对比图输出目录，默认 pipeline/figures")
    ap.add_argument("--no-pdf", action="store_true")
    a = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import mdpipe

    if a.compare:
        jobs = [mdpipe.job_dir(n) for n in a.compare]
        outdir = Path(a.outdir) if a.outdir else mdpipe.PIPE / "figures"
        plot_compare(jobs, a.compare, outdir, not a.no_pdf)
        return
    job = Path(a.dir) if a.dir else mdpipe.job_dir(a.job)
    if not job.exists():
        sys.exit(f"目录不存在：{job}")
    plot_job(job, a.name or job.name, not a.no_pdf)


if __name__ == "__main__":
    main()
