#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""画安装包行为流程图（纯 matplotlib，不依赖 graphviz）。

    python scripts/plot_installer_flow.py
    → pipeline/figures/installer_flow.png / .pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "pipeline" / "figures"

for fam in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "DejaVu Sans"):
    if any(fam in f.name for f in matplotlib.font_manager.fontManager.ttflist):
        matplotlib.rcParams["font.sans-serif"] = [fam]
        break
matplotlib.rcParams["axes.unicode_minus"] = False

INK = "#24303D"
MUTED = "#6B7785"
BLUE = "#DCEAF7"
BLUE_E = "#4A7FB5"
YELLOW = "#FFF3D6"
YELLOW_E = "#C9971E"
GREEN = "#E4F2E7"
GREEN_E = "#3E8E52"
RED = "#FBE6E3"
RED_E = "#B4483A"
GREY = "#F4F5F7"
GREY_E = "#AAB2BD"


def box(ax, cx, cy, w, h, text, fc=GREY, ec=GREY_E, fs=10.5, bold=False):
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.25,rounding_size=1.6",
        linewidth=1.4, edgecolor=ec, facecolor=fc, zorder=3))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs,
            color=INK, zorder=4,
            fontweight="bold" if bold else "normal", linespacing=1.35)


def diamond(ax, cx, cy, w, h, text, fs=9.5):
    ax.add_patch(Polygon([[cx, cy + h / 2], [cx + w / 2, cy],
                          [cx, cy - h / 2], [cx - w / 2, cy]],
                         closed=True, facecolor=YELLOW, edgecolor=YELLOW_E,
                         linewidth=1.4, zorder=3))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs,
            color=INK, zorder=4, linespacing=1.3)


def arrow(ax, p1, p2, label=None, rad=0.0, color=INK, ls="-", lw=1.35,
          lx=0.0, ly=0.0, fs=9):
    ax.add_patch(FancyArrowPatch(
        p1, p2, arrowstyle="-|>", mutation_scale=13, linewidth=lw,
        color=color, linestyle=ls, zorder=2,
        connectionstyle=f"arc3,rad={rad}", shrinkA=1, shrinkB=1))
    if label:
        ax.text((p1[0] + p2[0]) / 2 + lx, (p1[1] + p2[1]) / 2 + ly, label,
                ha="center", va="center", fontsize=fs, color=MUTED,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none"),
                zorder=5)


def poly_arrow(ax, pts, color=INK, ls="-"):
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color, lw=1.35,
            ls=ls, zorder=2, solid_capstyle="round")
    ax.add_patch(FancyArrowPatch(pts[-2], pts[-1], arrowstyle="-|>",
                                 mutation_scale=13, color=color, lw=1.35,
                                 linestyle=ls, zorder=2, shrinkA=0, shrinkB=1))


def main():
    fig, ax = plt.subplots(figsize=(11.2, 15.6))
    ax.set_xlim(-2, 120)
    ax.set_ylim(-24, 168)
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

    ax.text(60, 165, "ClawsGO MD 流水线 · 安装包行为",
            ha="center", va="center", fontsize=17, fontweight="bold", color=INK)
    ax.text(60, 161, "MDpipe-Setup-1.7.exe（Inno Setup 6）",
            ha="center", va="center", fontsize=10.5, color=MUTED)

    # 主流程
    box(ax, 30, 152, 46, 7.5, "双击 MDpipe-Setup-1.7.exe", BLUE, BLUE_E,
        fs=11.5, bold=True)
    box(ax, 30, 140, 46, 8, "选择安装位置\n默认 %LOCALAPPDATA%\\ClawsGO\\MDpipeline")
    diamond(ax, 30, 126, 42, 12, "检测到 Python 3.9+ ?")
    box(ax, 30, 109, 46, 9,
        "复制文件到安装目录\nscripts / mdp_templates / pipeline 模板\n"
        "tests / docs / 三个启动脚本")
    box(ax, 30, 96, 46, 8,
        "建运行目录\npipeline\\inbox | jobs | figures | splits、logs")
    box(ax, 30, 84, 46, 7.5,
        "按勾选创建\n桌面快捷方式 · 用户 PATH（可选）")
    box(ax, 30, 71, 46, 9,
        "写开始菜单\n图形界面 / 环境自检 / 命令行\n操作手册 / 环境清单 / 卸载")
    box(ax, 30, 57, 46, 6.5, "安装完成", GREEN, GREEN_E, bold=True)
    box(ax, 30, 45, 46, 7.5, "自动运行「环境自检」\n（只读 doctor，不装任何东西）",
        BLUE, BLUE_E)
    box(ax, 30, 32, 46, 8.5,
        "逐项检测 → 写 doctor_report.txt\n√ 可用    ! 可选缺失    × 必需缺失")
    diamond(ax, 30, 18, 40, 10.5, "有必需项 × ?")
    box(ax, 30, 3, 46, 7.5, "就绪：启动图形界面 / 投递作业", GREEN, GREEN_E,
        bold=True)

    # 缺 Python 分支
    box(ax, 86, 140, 46, 8.5,
        "弹窗：未检测到 Python\n（提醒装 3.9+ 并勾选 PATH）", RED, RED_E)
    diamond(ax, 86, 125, 32, 9.5, "仍要继续安装？", fs=9)
    box(ax, 86, 110, 34, 6.5, "退出安装程序", RED, RED_E)

    # 缺必需组件的补救
    box(ax, 86, 18, 42, 9.5,
        "按《安装与环境清单》补装\npip / winget / WSL / GROMACS …", YELLOW,
        YELLOW_E)

    # 箭头
    arrow(ax, (30, 148.2), (30, 144.2))
    arrow(ax, (30, 136), (30, 132.2))
    arrow(ax, (30, 120), (30, 113.7), "有", ly=1.6)
    arrow(ax, (51, 126), (63, 138.6), rad=-0.15, label="没有", lx=6, ly=3.5)
    arrow(ax, (86, 135.7), (86, 129.9))
    arrow(ax, (78, 121.2), (53, 110.5), rad=0.12, label="继续", lx=-4, ly=3.5)
    arrow(ax, (86, 120.2), (86, 113.4), "取消", lx=4.2, ly=0)
    arrow(ax, (30, 104.5), (30, 100.2))
    arrow(ax, (30, 92), (30, 87.9))
    arrow(ax, (30, 80.2), (30, 75.7))
    arrow(ax, (30, 66.5), (30, 60.4))
    arrow(ax, (30, 53.7), (30, 48.9))
    arrow(ax, (30, 41.2), (30, 36.5))
    arrow(ax, (30, 27.7), (30, 23.4))
    arrow(ax, (30, 12.7), (30, 6.9), "没有", lx=3.2, ly=0)
    arrow(ax, (50, 18), (64.8, 18), "有", ly=1.6)
    poly_arrow(ax, [(86, 22.9), (86, 45), (53.2, 45)])

    # 图例
    lx, ly = 8, -6
    for i, (label, fc, ec) in enumerate([
            ("开始 / 完成", GREEN, GREEN_E),
            ("安装动作", GREY, GREY_E),
            ("需要判断", YELLOW, YELLOW_E),
            ("异常 / 退出", RED, RED_E)]):
        ax.add_patch(FancyBboxPatch((lx + i * 27, ly - 1.6), 4.5, 3.2,
                                    boxstyle="round,pad=0.1,rounding_size=0.6",
                                    fc=fc, ec=ec, lw=1.2))
        ax.text(lx + i * 27 + 6, ly, label, va="center", fontsize=9.5,
                color=MUTED)

    # 安装包不做什么
    ax.add_patch(FancyBboxPatch((6, -22.5), 108, 11.5,
                                boxstyle="round,pad=0.3,rounding_size=1.2",
                                fc="#FBFBFC", ec=GREY_E, lw=1.2, ls=(0, (4, 3))))
    ax.text(9, -14.2,
            "安装包不做什么：不安装 Python / GROMACS / 显卡驱动 / WSL / conda；"
            "不下载任何运行时；\n"
            "不覆盖已存在的 config.json；卸载时保留 pipeline\\jobs 里的计算结果。",
            fontsize=9.5, color=MUTED, va="top", linespacing=1.6)

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"installer_flow.{ext}", dpi=220,
                    facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"已写出 {OUT / 'installer_flow.png'} 和 .pdf")


if __name__ == "__main__":
    sys.exit(main())
