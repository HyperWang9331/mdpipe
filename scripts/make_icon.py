#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""画 MDpipe 的应用图标（纯 matplotlib，不依赖设计软件）。

    python scripts/make_icon.py
    → pipeline/figures/appicon.png  (1024 px)
    → pipeline/figures/appicon.ico  (16–256 px 多尺寸，给 Windows 用)

图形取意：一条蛋白螺旋抱着一个小分子。螺旋用三维投影画，靠深浅表现前后；
右下角的琥珀色六边形是配体。缩到 16 px 时仍能认出“螺旋 + 小方块”。
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Polygon

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "pipeline" / "figures"

BG_TOP = "#122E4C"
BG_BOT = "#2B6C93"
HELIX_FRONT = "#F2FAFF"
HELIX_BACK = "#79A9C6"
AMBER = "#F5A93C"
AMBER_EDGE = "#FFD89A"


def draw(path: Path) -> None:
    fig = plt.figure(figsize=(4, 4), dpi=256)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")

    # 圆角底板 + 对角渐变
    tile = FancyBboxPatch((0.02, 0.02), 0.96, 0.96,
                          boxstyle="round,pad=0,rounding_size=0.19",
                          facecolor="none", edgecolor="none")
    ax.add_patch(tile)

    n = 512
    u = np.linspace(0, 1, n)
    grad = 0.55 * u[None, :] + 0.45 * (1 - u[:, None])
    cmap = LinearSegmentedColormap.from_list("md", [BG_TOP, BG_BOT])
    im = ax.imshow(grad, extent=(0, 1, 0, 1), origin="lower",
                   cmap=cmap, interpolation="bicubic", zorder=1)
    im.set_clip_path(tile)

    # 左上角柔光，避免大片死板深色（直接给 RGBA，免得 colormap 留下色带）
    yy, xx = np.mgrid[0:n, 0:n] / (n - 1)
    glow = np.exp(-(((xx - 0.24) ** 2 + (yy - 0.90) ** 2) / 0.09))
    rgba = np.empty((n, n, 4))
    rgba[..., :3] = 1.0
    rgba[..., 3] = 0.20 * glow ** 1.4
    ax.imshow(rgba, extent=(0, 1, 0, 1), origin="lower",
              interpolation="bicubic", zorder=2).set_clip_path(tile)

    # 螺旋：三维螺旋投影到二维，depth = cos 用颜色和线宽区分前后
    turns = 4.5
    t = np.linspace(0, 1, 900)
    x = 0.15 + 0.68 * t
    phase = 2 * np.pi * turns * t
    y = 0.535 + 0.128 * np.sin(phase)
    depth = np.cos(phase)                      # +1 在前，-1 在后

    pts = np.column_stack([x, y])
    segs = np.stack([pts[:-1], pts[1:]], axis=1)
    mid = 0.5 * (depth[:-1] + depth[1:])
    front = np.clip((mid + 1) / 2, 0, 1)

    def hexmix(a, b, k):
        ca = np.array([int(a[i:i + 2], 16) for i in (1, 3, 5)]) / 255
        cb = np.array([int(b[i:i + 2], 16) for i in (1, 3, 5)]) / 255
        return ca * k + cb * (1 - k)

    colors = np.array([hexmix(HELIX_FRONT, HELIX_BACK, k) for k in front])

    # 投影阴影
    ax.add_collection(LineCollection(
        segs + np.array([0.010, -0.012]), colors="black", alpha=0.18,
        linewidths=4.2 + 7.6 * front, capstyle="round", zorder=3))
    ax.add_collection(LineCollection(
        segs, colors=colors, linewidths=4.2 + 7.6 * front,
        capstyle="round", zorder=4))

    # 配体：琥珀色六边形，压在螺旋右端下方
    cx, cy, r = 0.760, 0.235, 0.098
    ang = np.deg2rad(np.arange(6) * 60 + 30)
    hexpts = np.column_stack([cx + r * np.cos(ang), cy + r * np.sin(ang)])
    ax.add_patch(Polygon(hexpts + np.array([0.008, -0.010]), closed=True,
                         facecolor="black", alpha=0.18, edgecolor="none",
                         zorder=3))
    ax.add_patch(Polygon(hexpts, closed=True, facecolor=AMBER,
                         edgecolor=AMBER_EDGE, linewidth=2.4, zorder=5))
    # 六边形内留白，像简化结构式
    ax.add_patch(Polygon(hexpts * 0.55 + np.array([cx, cy]) * 0.45, closed=True,
                         facecolor="#FFF6E6", edgecolor="none", zorder=6))

    fig.savefig(path, dpi=256, transparent=True)
    plt.close(fig)


def make_ico(png: Path, ico: Path) -> None:
    from PIL import Image
    im = Image.open(png).convert("RGBA")
    sizes = [(s, s) for s in (16, 24, 32, 48, 64, 128, 256)]
    im.save(ico, format="ICO", sizes=sizes)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "appicon.png"
    ico = OUT / "appicon.ico"
    draw(png)
    make_ico(png, ico)
    print(f"已写出 {png} 和 {ico}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
