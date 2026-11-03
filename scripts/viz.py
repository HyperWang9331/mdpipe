#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""viz.py —— 复合物可视化与三条补充轨迹指标。

补的是六联图里没有、但审稿人常问的东西：

1. **复合物 RMSD**：Protein_LIG 整体的 RMSD，说明"蛋白和配体作为一个复合物"
   有没有整体漂移（骨架 RMSD 只说明蛋白自己）。
2. **埋藏 SASA**：SASA(蛋白单独) − SASA(蛋白在复合物里)，即配体盖住的那部分
   表面积。结合得越紧，这条曲线越稳定、越大。
3. **质心距离**：蛋白质心到配体质心的距离，配体离开口袋会很直观地抬升。

再加两张结构图（PyMOL 渲染）：
- 0 / 中 / 末三个时刻的快照，三张图共用同一个视角，可以直接比较；
- 逐残基 RMSF 映射到结构上的 B-factor 着色图。

除了 PNG，同时写出 `.pml` 脚本和 `.pse` 会话文件，方便手工改视角、换配色、
加标注后重新渲染。

单独用法
--------
    python scripts/viz.py --job Gal_1
    python scripts/viz.py --dir path/to/job --frames 0 25 50 75 100
    python scripts/viz.py --job Gal_1 --no-pymol          # 只算三条曲线
"""
from __future__ import annotations

import argparse
import contextlib
import io
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib import rcParams                       # noqa: E402
from matplotlib.ticker import FuncFormatter            # noqa: E402

from make_figures import load_xvg, save, style_panel, RED, BLUE  # noqa: E402

rcParams["font.family"] = ["Arial", "Microsoft YaHei", "DejaVu Sans"]
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 9.5
rcParams["axes.linewidth"] = 0.9
rcParams["axes.grid"] = False

PALETTE = [RED, BLUE, "#008000", "#FF7F00"]
DPI_PNG = 400


# --------------------------------------------------------------------------
# gmx 调用
# --------------------------------------------------------------------------
def _gmx(cfg, args, cwd: Path, stdin_text: str | None = None, log: Path | None = None):
    cmd = [cfg.get("gmx", "gmx")] + [str(a) for a in args]
    p = subprocess.run(cmd, cwd=str(cwd), input=stdin_text, text=True,
                       capture_output=True, errors="replace")
    if log is not None:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text((p.stdout or "") + "\n" + (p.stderr or ""), encoding="utf-8")
    if p.returncode != 0:
        tail = ((p.stderr or "") + (p.stdout or ""))[-900:]
        raise RuntimeError(f"gmx {' '.join(map(str, args))} 失败：\n{tail}")
    return p.stdout


def _traj_ns(ana: Path) -> float:
    for pat in ("rmsd_backbone_*.xvg", "rmsd_backbone.xvg"):
        for p in sorted(ana.glob(pat)):
            x, _ = load_xvg(p, 1)
            if x is not None and len(x):
                return float(x[-1])
    return 0.0


def ensure_inputs(cfg, job: Path, meta: dict, rep: int) -> dict:
    """找到（必要时生成）viz 需要的轨迹和参考结构。

    分析目录里的 md_prot.xtc 是全帧的，直接用；如果它被归档删掉了，就在
    viz/traj/ 下按 100 ps 抽一版——三条指标对 10 帧/ns 的采样完全够用，
    体积只有原轨迹的 1/5。
    """
    d = job
    top = d / "02_topology"
    sfx = f"_r{rep}" if meta.get("replicates", 1) > 1 else ""
    prod = d / f"06_prod{sfx}"
    ana = d / f"07_analysis{sfx}"
    out = d / "08_figures" / "viz"
    out.mkdir(parents=True, exist_ok=True)
    tpr, xtc, ndx = prod / "md.tpr", prod / "md.xtc", top / "index.ndx"
    if not tpr.exists():
        raise RuntimeError(f"缺 {tpr}")
    res = meta.get("resname", "LIG")

    prot = ana / "md_prot.xtc"
    ref = ana / "ref0.gro"
    if not (prot.exists() and ref.exists()):
        tdir = out / "traj"
        tdir.mkdir(parents=True, exist_ok=True)
        prot = tdir / f"md_prot{sfx}.xtc"
        ref = tdir / f"ref0{sfx}.gro"
        if not prot.exists():
            _gmx(cfg, ["trjconv", "-s", tpr, "-f", xtc, "-n", ndx,
                       "-o", prot, "-pbc", "mol", "-center", "-dt", "100"],
                 cwd=prod, stdin_text="Protein_LIG\nProtein_LIG\n",
                 log=out / f"_trjconv_prot{sfx}.log")
        if not ref.exists():
            _gmx(cfg, ["trjconv", "-s", tpr, "-f", prot, "-n", ndx,
                       "-o", ref, "-dump", "0"],
                 cwd=prod, stdin_text="Protein_LIG\n",
                 log=out / f"_ref0{sfx}.log")
    return {"job": d, "ana": ana, "prod": prod, "tpr": tpr, "xtc": xtc,
            "rep": rep,
            "ndx": ndx, "ref": ref, "prot": prot, "out": out, "res": res,
            "sfx": sfx, "ns": _traj_ns(ana)}


# --------------------------------------------------------------------------
# 三条补充指标
# --------------------------------------------------------------------------
def metric_complex_rmsd(cfg, m: dict) -> Path:
    out = m["out"] / f"complex_rmsd{m['sfx']}.xvg"
    _gmx(cfg, ["rms", "-s", m["ref"], "-f", m["prot"], "-n", m["ndx"],
               "-o", out, "-tu", "ns"],
         cwd=m["prod"], stdin_text="Protein_LIG\nProtein_LIG\n",
         log=m["out"] / f"_complex_rmsd{m['sfx']}.log")
    return out


def metric_com_distance(cfg, m: dict) -> Path:
    out = m["out"] / f"com_distance{m['sfx']}.xvg"
    _gmx(cfg, ["distance", "-s", m["ref"], "-f", m["prot"], "-n", m["ndx"],
               "-select", 'com of group "Protein" plus com of group "LIG"',
               "-oall", out, "-tu", "ns"],
         cwd=m["prod"], log=m["out"] / f"_com_distance{m['sfx']}.log")
    return out


def metric_buried_sasa(cfg, m: dict) -> dict:
    """两个埋藏指标：

    * 埋藏面积 = SASA(蛋白单独) − SASA(蛋白在复合物里)，单位 nm²；
    * 埋藏比例 = 1 − SASA(配体在复合物里) / SASA(配体单独)，即配体被包住的
      表面积百分比。文献里常用后者，>85% 一般认为配体埋得足够好。

    SASA 对帧不敏感，统一按 0.2 ns 抽帧算（100 ns 只有 500 帧），比逐帧快 10 倍。
    """
    sfx, res = m["sfx"], m["res"]
    base = ["sasa", "-s", m["ref"], "-f", m["prot"], "-n", m["ndx"],
            "-tu", "ns", "-dt", "0.2"]
    # 复合物表面，分解出蛋白和配体各自的面积（两次运行，-output 一次只能给一组）
    run_p = m["out"] / f"_sasa_incomplex_p{sfx}.xvg"
    _gmx(cfg, base + ["-surface", f'group "Protein" or group "{res}"',
                      "-output", 'group "Protein"', "-o", run_p],
         cwd=m["prod"], log=m["out"] / f"_sasa_p{sfx}.log")
    run_l = m["out"] / f"_sasa_incomplex_l{sfx}.xvg"
    _gmx(cfg, base + ["-surface", f'group "Protein" or group "{res}"',
                      "-output", f'group "{res}"', "-o", run_l],
         cwd=m["prod"], log=m["out"] / f"_sasa_l{sfx}.log")
    # 配体单独的面积
    run_la = m["out"] / f"_sasa_lig_alone{sfx}.xvg"
    _gmx(cfg, base + ["-surface", f'group "{res}"', "-o", run_la],
         cwd=m["prod"], log=m["out"] / f"_sasa_la{sfx}.log")

    t_c, (col_prot,) = _read_xvg_cols(run_p, (2,))
    _, (col_lig,) = _read_xvg_cols(run_l, (2,))
    _, (col_lig_alone,) = _read_xvg_cols(run_la, (1,))

    alone = m["ana"] / f"sasa_{m['rep']}.xvg"
    if not alone.exists():
        cands = sorted(m["ana"].glob("sasa_*.xvg"))
        alone = cands[0] if cands else None
    t_a, col_alone = None, None
    if alone is not None and alone.exists():
        t_a, (col_alone,) = _read_xvg_cols(alone, (1,))
    if col_alone is None:                    # 老作业没有单独面积就现场算一遍
        tmp = m["out"] / f"_sasa_alone{sfx}.xvg"
        _gmx(cfg, base + ["-o", tmp], cwd=m["prod"], stdin_text="Protein\n",
             log=m["out"] / f"_sasa_alone{sfx}.log")
        t_a, (col_alone,) = _read_xvg_cols(tmp, (1,))

    if t_a is not None and (len(t_a) != len(t_c) or not np.allclose(t_a, t_c)):
        col_alone = np.interp(t_c, t_a, col_alone)      # 时间轴不一致就插值
    buried = col_alone - col_prot
    frac = 100.0 * (1.0 - col_lig / np.maximum(col_lig_alone, 1e-9))

    out1 = m["out"] / f"buried_sasa{sfx}.xvg"
    _write_xy(out1, t_c, buried, "Buried SASA", "Area (nm\\S2\\N)",
              "# 埋藏面积 = SASA(蛋白单独) - SASA(蛋白在复合物中)")
    out2 = m["out"] / f"ligand_buried{sfx}.xvg"
    _write_xy(out2, t_c, frac, "Ligand buried fraction", "Buried (%)",
              "# 埋藏比例 = 1 - SASA(配体在复合物中)/SASA(配体单独)")
    return {"buried_sasa": out1, "ligand_buried": out2}


def _write_xy(path: Path, t, y, title, ylab, comment):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(comment + "\n# created by viz.py\n")
        fh.write(f'@    title "{title}"\n@    xaxis  label "Time (ns)"\n'
                 f'@    yaxis  label "{ylab}"\n@TYPE xy\n')
        for ti, yi in zip(t, y):
            fh.write(f"{ti:.4f} {yi:.4f}\n")


def _read_xvg_cols(path: Path, cols):
    """读 xvg 的若干列（cols 从 1 开始，0 是时间）。"""
    t, vals = [], [[] for _ in cols]
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s or s[0] in "#@":
            continue
        q = s.split()
        try:
            t.append(float(q[0]))
            for j, c in enumerate(cols):
                vals[j].append(float(q[c]))
        except (ValueError, IndexError):
            continue
    return np.array(t), [np.array(v) for v in vals]


def compute_metrics(cfg, m: dict) -> dict:
    print(f"  复合物 RMSD …", flush=True)
    complex_rmsd = metric_complex_rmsd(cfg, m)
    print(f"  质心距离 …", flush=True)
    comdist = metric_com_distance(cfg, m)
    print(f"  埋藏 SASA …（这条最慢）", flush=True)
    buried = metric_buried_sasa(cfg, m)
    return {"complex_rmsd": complex_rmsd, "com_distance": comdist, **buried}


def plot_metrics(curves: dict, outbase: Path, name: str, want_pdf=True):
    """四面板：复合物 RMSD / 埋藏面积 / 配体埋藏比例 / 质心距离。"""
    panels = [("Complex RMSD", "complex_rmsd", "Complex RMSD (nm)"),
              ("Buried SASA", "buried_sasa", "Buried area (nm$^2$)"),
              ("Ligand burial", "ligand_buried", "Buried ligand SASA (%)"),
              ("COM distance", "com_distance", "Protein–ligand COM (nm)")]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.2))
    stats = []
    for i, (title, key, ylab) in enumerate(panels):
        ax = axes[i // 2][i % 2]
        series = curves.get(key) or []
        for k, (x, y) in enumerate(series):
            ax.plot(x, y, color=PALETTE[k % len(PALETTE)], lw=0.6,
                    label=f"Replicate {k + 1}" if len(series) > 1 else None)
            tail = y[int(len(y) * 0.5):]
            stats.append((title, k + 1, float(y.mean()), float(y.std()),
                          float(tail.mean()), float(tail.std())))
        if key == "ligand_buried":
            ax.axhline(85, color="#4a5568", lw=0.7, ls="--", zorder=0)
            ax.text(0.985, 85, "85%", transform=ax.get_yaxis_transform(),
                    ha="right", va="bottom", fontsize=7, color="#4a5568")
            ax.set_ylim(0, 100)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel(ylab)
        style_panel(ax, "ABCD"[i])
        if len(series) > 1:
            ax.legend(frameon=False, fontsize=7.5, loc="best",
                      handlelength=1.6, borderaxespad=0.3)
    plt.subplots_adjust(wspace=0.30, hspace=0.42)
    outbase.parent.mkdir(parents=True, exist_ok=True)
    outs = save(fig, outbase, want_pdf)
    plt.close(fig)
    print("  " + ", ".join(p.name for p in outs))
    return stats


# --------------------------------------------------------------------------
# PyMOL 场景
# --------------------------------------------------------------------------
def _pymol():
    import pymol
    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd
    cmd.feedback("disable", "all")
    return cmd


@contextlib.contextmanager
def _silent():
    """同时压住 PyMOL 的 C 级回显和 Python 打印。"""
    buf = io.StringIO()
    try:
        sys.stdout.flush()
        sys.stderr.flush()
        fd_out, fd_err = os.dup(1), os.dup(2)
        null = os.open(os.devnull, os.O_WRONLY)
        os.dup2(null, 1)
        os.dup2(null, 2)
    except Exception:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            yield
        return
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            yield
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(fd_out, 1)
        os.dup2(fd_err, 2)
        os.close(null)
        os.close(fd_out)
        os.close(fd_err)


def snapshot_lines(pdbs, times, res: str) -> list[str]:
    """一帧一个 object，全部叠在同一坐标系里；只切换显示，不重设视角。"""
    lines = ["# 复合物快照：半透明卡通 + 配体 sticks + 口袋侧链 sticks",
             "# 样式参考 PyMOL 社区常用的论文图设置：白底 + ray_trace_mode 1 + 卡通透明度",
             "# 每帧一个 object（t0ns / t50ns / …），用 enable/disable 切换"]
    for t, pdb in zip(times, pdbs):
        lines.append(f'load "{pdb.as_posix()}", t{t:g}ns')
    lines += [
        "hide everything",
        "show cartoon, polymer",
        "color grey70, polymer",
        "set cartoon_transparency, 0.60",
        "set cartoon_side_chain_helper, on",
        "hide (hydro and polymer)",
        f'select pocket, byres (polymer within 4.0 of resn {res}) and not (name N+C+O+CA)',
        "show sticks, pocket",
        "color grey50, pocket and elem C",
        "util.cnc pocket",
        f'show sticks, resn {res}',
        f'color yellow, resn {res} and elem C',
        f'util.cnc resn {res}',
        "set stick_radius, 0.16, pocket",
        f'set stick_radius, 0.24, resn {res}',
        "set ray_trace_mode, 1",
        "set ray_trace_gain, 0",
        "set ray_shadows, 0",
        "set antialias, 2",
        "set specular, 0.15",
        "set depth_cue, 0",
        "bg_color white",
        "disable all",
        f"enable t{times[0]:g}ns",
        f"orient resn {res} and t{times[0]:g}ns",
        f"zoom resn {res} and t{times[0]:g}ns, 10",
    ]
    return lines


def render_snapshots(cmd, pdbs, times, res, outdir: Path, prefix: str,
                     pse_name: str, width=760, height=660):
    """一张会话、固定视角、逐帧切 object，保证三张图可以直接比较。"""
    lines = snapshot_lines(pdbs, times, res)
    pngs = []
    with _silent():
        cmd.reinitialize()
        for ln in lines:
            cmd.do(ln)
        for t in times:
            cmd.do("disable all")
            cmd.do(f"enable t{t:g}ns")
            # 只重设中心与距离，旋转矩阵沿用首帧 → 朝向一致；取景包含配体
            # 和它周围 12 Å 的蛋白，配体跑出口袋时也不会只剩一个孤零零的分子
            cmd.do(f"zoom (byres (polymer within 12 of resn {res}) or resn {res}) "
                   f"and t{t:g}ns, 2")
            png = outdir / f"_scene{prefix}_{t:g}ns.png"
            cmd.png(str(png), width=width, height=height, dpi=150, ray=1)
            pngs.append(png)
        cmd.do("disable all")
        cmd.do(f"enable t{times[0]:g}ns")
        cmd.save(str(outdir / pse_name))
    return pngs, lines


def scene_bfactor(pdb: Path, ligand_pdb: Path | None, res: str,
                  blo: float, bhi: float, out_png: Path, width=760, height=660):
    cmd = scene_bfactor.cmd
    cmd.reinitialize()
    lines = [
        "# B-factor 着色：颜色 = 逐残基 RMSF",
        f'load "{pdb.as_posix()}", bfac',
        "hide everything",
        "show cartoon, polymer",
        f"spectrum b, blue_white_red, minimum={blo:.1f}, maximum={bhi:.1f}",
        "set cartoon_transparency, 0.0",
        "set cartoon_side_chain_helper, on",
    ]
    if ligand_pdb is not None:
        lines += [f'load "{ligand_pdb.as_posix()}", lig',
                  f'hide everything, lig',
                  f'show sticks, lig and resn {res}',
                  f'color yellow, lig and resn {res} and elem C',
                  f'util.cnc lig and resn {res}']
    lines += [
        "set stick_radius, 0.2",
        "set ray_trace_mode, 0",       # 着色图不要黑边，颜色才是主角
        "set ray_trace_gain, 0",
        "set ray_shadows, 0",
        "set antialias, 2",
        "set specular, 0.05",
        "set depth_cue, 0",
        "bg_color white",
        "orient polymer",
    ]
    for ln in lines:
        cmd.do(ln)
    cmd.png(str(out_png), width=width, height=height, dpi=150, ray=1)
    return out_png, lines


def _write_pml(path: Path, lines: list[str], pse: Path | None = None,
               title: str = ""):
    body = (["# " + title, "# 这个脚本由 viz.py 自动生成，可以直接编辑后重跑：",
             f"#   pymol {path.name}", ""] + lines)
    if pse is not None:
        body += [f'# 会话文件：{pse.name}（在 PyMOL 里 File → Open 打开）']
    path.write_text("\n".join(body) + "\n", encoding="utf-8")


def _compose_snapshots(pngs, times, outbase: Path, want_pdf=True):
    fig, axes = plt.subplots(1, len(pngs), figsize=(4.3 * len(pngs), 4.4))
    if len(pngs) == 1:
        axes = [axes]
    for ax, png, t in zip(axes, pngs, times):
        ax.imshow(plt.imread(png))
        ax.axis("off")
        ax.set_title(f"t = {t:g} ns", fontsize=11, pad=6)
    plt.subplots_adjust(wspace=0.02)
    outs = [outbase.with_suffix(".png")]
    if want_pdf:
        outs += [outbase.with_suffix(".pdf")]
    for p in outs:
        fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white",
                    pad_inches=0.08)
    plt.close(fig)
    print("  " + ", ".join(p.name for p in outs))
    return outs


def _compose_bfactor(png: Path, blo: float, bhi: float, outbase: Path, name: str,
                     want_pdf=True):
    import matplotlib.cm as cm
    from matplotlib.colors import Normalize
    fig = plt.figure(figsize=(5.6, 4.6))
    ax = fig.add_axes([0.02, 0.06, 0.78, 0.9])
    ax.imshow(plt.imread(png))
    ax.axis("off")
    cax = fig.add_axes([0.83, 0.12, 0.035, 0.76])
    sm = cm.ScalarMappable(norm=Normalize(blo, bhi), cmap="coolwarm")
    cb = fig.colorbar(sm, cax=cax, format=FuncFormatter(lambda v, _: f"{v / 1000:g}"))
    cb.set_label("RMSF (nm)", fontsize=10)
    cb.ax.tick_params(labelsize=8)
    short = name.split("__")[-1].replace("_", " ")
    fig.text(0.02, 0.965, f"{short} — RMSF mapped on structure", fontsize=10)
    fig.text(0.02, 0.935, f"colour scale capped at the 90th percentile "
                          f"({bhi / 1000:.2f} nm); higher values saturate red",
             fontsize=7.5, color="#4a5568")
    outs = [outbase.with_suffix(".png")]
    if want_pdf:
        outs += [outbase.with_suffix(".pdf")]
    for p in outs:
        fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  " + ", ".join(p.name for p in outs))
    return outs


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def run_viz(cfg, job: Path, name: str, meta: dict, reps=None,
            frames=None, want_pdf=True, do_pymol=True):
    job = Path(job)
    nrep = int(meta.get("replicates", 1))
    reps = reps or list(range(1, nrep + 1))
    out = job / "08_figures" / "viz"
    out.mkdir(parents=True, exist_ok=True)

    curves = {"complex_rmsd": [], "buried_sasa": [], "ligand_buried": [],
              "com_distance": []}
    frames_by_rep = {}
    bfac = {}
    for rep in reps:
        sfx = f"_r{rep}" if nrep > 1 else ""
        print(f"[{name}] 重复 {rep}：", flush=True)
        try:
            m = ensure_inputs(cfg, job, meta, rep)
        except Exception as exc:
            print(f"  跳过：{exc}")
            continue
        paths = compute_metrics(cfg, m)
        for key, p in paths.items():
            x, y = load_xvg(p, 1)
            if x is not None:
                curves[key].append((x, y))

        # 快照时间：默认 0 / 中 / 末（轨迹不到 100 ns 就按实际长度取）
        total = m["ns"] or (_traj_ns(m["ana"]) if m["ana"].exists() else 0) or 100.0
        times = frames or [0.0, round(total / 2, 1), round(total, 1)]
        times = [float(t) for t in times]
        pdbs = []
        for t in times:
            fp = out / f"frames{sfx}" / f"frame_{t:g}ns.pdb"
            fp.parent.mkdir(parents=True, exist_ok=True)
            if not fp.exists():
                # -fit rot+trans：把每一帧按骨架叠合到首帧，否则蛋白在盒子里
                # 自由转动，同一个视角下三个时刻会朝向完全不同，没法比较。
                _gmx(cfg, ["trjconv", "-s", m["ref"], "-f", m["prot"], "-n", m["ndx"],
                           "-o", fp, "-dump", f"{t * 1000:.0f}", "-fit", "rot+trans"],
                     cwd=m["prod"], stdin_text="Backbone\nProtein_LIG\n",
                     log=out / f"_frame_{t:g}ns{sfx}.log")
            pdbs.append(fp)
        frames_by_rep[rep] = (times, pdbs, m)

        # B-factor：gmx rmsf -oq 直接把 RMSF 写进 B 列
        bp = out / f"bfactor{sfx}.pdb"
        if not bp.exists():
            _gmx(cfg, ["rmsf", "-s", m["ref"], "-f", m["prot"], "-n", m["ndx"],
                       "-o", out / f"rmsf{sfx}.xvg", "-oq", bp],
                 cwd=m["prod"], stdin_text="Protein\n",
                 log=out / f"_rmsf{sfx}.log")
        bfac[rep] = (bp, m)

    if any(curves.values()):
        print(f"[{name}] 三条补充曲线 {out.name}/", flush=True)
        plot_metrics(curves, out / f"{name}_complex_metrics", name, want_pdf)

    if not do_pymol or not frames_by_rep:
        return

    try:
        with _silent():
            cmd = _pymol()
    except Exception as exc:                    # PyMOL 不在就只交曲线
        print(f"[{name}] PyMOL 不可用，跳过结构图：{exc}")
        return
    scene_bfactor.cmd = cmd
    res = meta.get("resname", "LIG")

    for rep, (times, pdbs, m) in frames_by_rep.items():
        sfx = f"_r{rep}" if nrep > 1 else ""
        print(f"[{name}] 重复 {rep}：渲染快照 {times}", flush=True)
        pngs, lines = render_snapshots(cmd, pdbs, times, res, out, f"{sfx}_",
                                       f"{name}_snapshots{sfx}.pse",
                                       width=760, height=660)
        _compose_snapshots(pngs, times, out / f"{name}_snapshots{sfx}", want_pdf)
        _write_pml(out / f"snapshots{sfx}.pml", lines,
                   title=f"{name} 重复{rep}：{times} ns 快照（enable/disable 切换）")
        print(f"  {out.name}/snapshots{sfx}.pml  +  {name}_snapshots{sfx}.pse")

        # B-factor 图
        bp, mm = bfac[rep]
        blo, bhi = _pdb_bstats(bp)
        if bhi <= 0:
            print("  B-factor 全为 0，跳过着色图")
            continue
        print(f"[{name}] 重复 {rep}：B-factor 着色（RMSF "
              f"0–{bhi / 1000:.2f} nm，90 分位封顶）", flush=True)
        lig = pdbs[0] if pdbs else None
        png = out / f"_scene_bfactor{sfx}.png"
        with _silent():
            _, lines = scene_bfactor(bp, lig, res, blo, bhi, png)
        _compose_bfactor(png, blo, bhi, out / f"{name}_bfactor{sfx}", name, want_pdf)
        _write_pml(out / f"bfactor{sfx}.pml", lines,
                   title=f"{name} 重复{rep}：RMSF → B-factor")
        try:
            with _silent():
                cmd.reinitialize()
                for ln in lines:
                    cmd.do(ln)
                cmd.save(str(out / f"{name}_bfactor{sfx}.pse"))
        except Exception as exc:
            print(f"  B-factor 会话没写成：{exc}")


def _pdb_bstats(path: Path):
    """gmx rmsf -oq 把 RMSF 以「nm × 1000」写进 B 列。色标用 0–90 分位：
    末端和 loop 上的少数原子动辄 0.3–0.8 nm，按最大值定标会把整个桶压成深蓝。"""
    vals = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                continue
    if not vals:
        return 0.0, 0.0
    return 0.0, float(np.percentile(vals, 90))


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description="复合物可视化（补充指标 + 结构图）")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--job", help="队列里的作业名（可写多个，用逗号分隔）")
    g.add_argument("--dir", help="任意作业目录")
    g.add_argument("--all", action="store_true")
    ap.add_argument("--replicates", nargs="+", type=int)
    ap.add_argument("--frames", nargs="+", type=float,
                    help="快照时间 (ns)，默认 0 / 中 / 末")
    ap.add_argument("--no-pymol", action="store_true", help="只算三条曲线")
    ap.add_argument("--no-pdf", action="store_true")
    a = ap.parse_args(argv)

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import mdpipe
    cfg = mdpipe.config()

    if a.dir:
        job = Path(a.dir)
        meta = mdpipe.load_json(job / "job.json", {})
        run_viz(cfg, job, meta.get("name") or job.name, meta,
                a.replicates, a.frames, not a.no_pdf, not a.no_pymol)
        return
    names = (mdpipe.queue()["jobs"] if a.all
             else [s.strip() for s in a.job.split(",") if s.strip()])
    for n in names:
        job = mdpipe.job_dir(n)
        meta = mdpipe.job_meta(n)
        if not job.exists() or not meta:
            print(f"跳过 {n}：没有作业目录或 job.json")
            continue
        run_viz(cfg, job, n, meta, a.replicates, a.frames,
                not a.no_pdf, not a.no_pymol)


if __name__ == "__main__":
    main()
