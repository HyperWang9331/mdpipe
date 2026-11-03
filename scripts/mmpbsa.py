#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""mmpbsa.py —— 批量 MM-PBSA 与能量分解，出 ΔG 对比图。

底层还是 gmx_MMPBSA（WSL 里的 conda 环境，路径由 config.json 指定）。包装脚本处理三件麻烦事：

1. 那个环境里的 gmx 是 2019 版，读不了 2026 生成的 tpr —— 先用它重建结构 tpr；
2. ParmEd 解析带 `#include` 的 topol.top 会把键参数读坏 —— 先用 2026 的
   `grompp -pp` 展开成平铺拓扑；
3. ParmEd 找不到 amber14sb.ff —— 必须给 GMXLIB；
4. 原始 md.xtc 里蛋白跨周期边界，cpptraj 会把它看成断裂分子，BOND 项飙到 1e7
   —— 先 `trjconv -pbc mol -center` 居中。

取样窗口默认取轨迹最后 20 ns（`--last-ns`，0 表示全程），与文献常用做法一致。

用法：
    python scripts/mmpbsa.py --job Gal_1 --replicates 1 2 --frames 10
    python scripts/mmpbsa.py --job Gal_1 --last-ns 50        # 换取样窗口
    python scripts/mmpbsa.py --job Gal_1 Glc_1 --plot        # 只重画已有结果

产物在作业的 `08_mmpbsa[_rN]/`；多体系对比图在 `pipeline/figures/`。
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

def _pipe_config() -> dict:
    """读 pipeline/config.json；路径统一在那里改，脚本里不再写死。"""
    p = Path(__file__).resolve().parent.parent / "pipeline" / "config.json"
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {}


def env_paths() -> dict:
    """WSL 侧的路径与发行版名，全部可由 config.json 覆盖。"""
    c = _pipe_config()
    gmx = c.get("gmx") or "gmx"
    share = c.get("gmx_share") or win_to_wsl(
        Path(gmx).resolve().parent.parent / "share" / "gromacs" / "top")
    return {
        "env": c.get("mmpbsa_env") or "",
        "gmxlib": share,
        "distro": c.get("wsl_distro") or "",
        "home": c.get("wsl_home") or "",
    }

# ΔTOTAL 之后这些项值得单独看
TERMS = ["VDWAALS", "EEL", "EPB", "ENPOLAR", "EDISPER", "GGAS", "GSOLV", "TOTAL"]


def win_to_wsl(p) -> str:
    s = str(p).replace("\\", "/")
    if len(s) > 1 and s[1] == ":":
        s = "/mnt/" + s[0].lower() + s[2:]
    return s


def _wsl(inner: str):
    e = env_paths()
    return ["wsl.exe", "-d", e["distro"], "-e", "env", "-i",
            "PATH=/usr/local/bin:/usr/bin:/bin", f"HOME={e['home']}",
            "bash", "--noprofile", "--norc", "-c", inner]


def _run(cmd, cwd=None, log=None, check=True):
    fh = open(log, "w", encoding="utf-8", errors="replace") if log else subprocess.DEVNULL
    try:
        r = subprocess.run(cmd, cwd=str(cwd) if cwd else None,
                           stdout=fh, stderr=subprocess.STDOUT, text=True)
    finally:
        if log:
            fh.close()
    if check and r.returncode != 0:
        raise RuntimeError(f"命令失败：{' '.join(map(str, cmd))}\n日志：{log}")
    return r.returncode


def ndx_groups(ndx: Path) -> list[str]:
    names = []
    for line in ndx.read_text(encoding="utf-8", errors="replace").splitlines():
        m = re.match(r"\[\s*(.+?)\s*\]", line)
        if m:
            names.append(m.group(1))
    return names


def group_index(ndx: Path, name: str) -> int:
    """gmx_MMPBSA 的 -cg 用 0-based 组号（和 gmx 的组选择提示一致）。"""
    names = ndx_groups(ndx)
    if name not in names:
        raise RuntimeError(f"{ndx.name} 里没有组 {name}（有：{', '.join(names[:12])}…）")
    return names.index(name)


def prepare_topology(job: Path, out: Path, gmx: str):
    top = job / "02_topology"
    mdp = out / "tpr2019.mdp"
    mdp.write_text("""integrator      = md
dt              = 0.002
nsteps          = 0
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
tcoupl          = no
pcoupl          = no
constraints     = none
gen-vel         = no
""", encoding="utf-8")
    _run([gmx, "grompp", "-f", str(mdp), "-c", str(top / "system.gro"),
          "-p", str(top / "topol.top"), "-n", str(top / "index.ndx"),
          "-pp", str(out / "processed.top"), "-o", str(out / "pp.tpr"),
          "-maxwarn", "2"], cwd=top, log=out / "grompp_pp.log")
    # 2019 的 gmx 重建 tpr
    e = env_paths()
    inner = (f'export PATH={e["env"]}/bin:$PATH; export GMXLIB={e["gmxlib"]}; '
             f'cd "{win_to_wsl(job)}" && gmx grompp -f {win_to_wsl(mdp)} '
             f'-c 02_topology/system.gro -p 02_topology/topol.top '
             f'-n 02_topology/index.ndx -o {win_to_wsl(out)}/complex.tpr -maxwarn 2')
    _run(_wsl(inner), log=out / "grompp2019.log")


def write_input(out: Path, frames: int, nframe: int, last_ps: float,
                last_ns: float, istrng: float, decomp: bool):
    # 取样窗口默认取轨迹最后 20 ns（平衡段），这是文献里最常见的做法
    dt = last_ps / max(1, nframe - 1)
    if last_ns and last_ns * 1000 < last_ps:
        start = max(1, int(round(nframe - last_ns * 1000 / dt)))
    else:
        start = 1
    span = nframe - start + 1
    interval = max(1, span // max(1, frames))
    text = f"""&general
  startframe={start}, endframe={nframe}, interval={interval},
  verbose=1, keep_files=0,
/
&pb
  istrng={istrng}, fillratio=4.0, radiopt=0,
/
"""
    if decomp:
        # 分解用 GB：PB 的逐残基分解要重算每个残基的溶剂化，慢一个数量级
        text += """&gb
  igb=5, saltcon=0.150,
/
&decomp
  idecomp=1, dec_verbose=1,
/
"""
    (out / "mmpbsa.in").write_text(text, encoding="utf-8")


def run_one(job: Path, rep: int | None, frames=10, istrng=0.15,
            decomp=False, gmx="gmx",
            resname="LIG", last_ns=20.0, keep_center=False) -> dict:
    job = Path(job).resolve()
    nrep_dirs = sorted(job.glob("06_prod*"))
    if rep is None:
        if (job / "06_prod").exists():
            prod, sfx = job / "06_prod", ""
        elif nrep_dirs:
            d = nrep_dirs[0]
            prod, sfx = d, d.name[len("06_prod"):]
        else:
            raise RuntimeError(f"{job} 下没有 06_prod*")
    else:
        sfx = f"_r{rep}" if len(nrep_dirs) > 1 or (job / f"06_prod_r{rep}").exists() else ""
        prod = job / f"06_prod{sfx}"
    if not (prod / "md.xtc").exists():
        raise RuntimeError(f"缺 {prod}/md.xtc")

    out = job / f"08_mmpbsa{sfx}"
    out.mkdir(exist_ok=True)
    ndx = job / "02_topology" / "index.ndx"

    whole = _whole_traj(gmx, prod, ndx, out, out / "run.log")
    nframe, last_ps = _n_frames(whole, gmx)
    prepare_topology(job, out, gmx)
    write_input(out, frames, nframe, last_ps, last_ns, istrng, decomp)

    prot = group_index(ndx, "Protein")
    lig = group_index(ndx, resname)
    e = env_paths()
    cmd = (f'export PATH={e["env"]}/bin:$PATH; export GMXLIB={e["gmxlib"]}; '
           f'cd "{win_to_wsl(out)}" && gmx_MMPBSA -O -i mmpbsa.in '
           f'-cs complex.tpr -ci {win_to_wsl(ndx)} -cg {prot} {lig} '
           f'-ct {win_to_wsl(whole)} -cp processed.top '
           f'-o FINAL_RESULTS.dat')
    _run(_wsl(cmd), log=out / "run.log")

    res = parse_results(out / "FINAL_RESULTS.dat")
    res.update({"job": job.name, "rep": sfx.lstrip("_r") or "1",
                "dir": str(out), "decomp": decomp,
                "frames": frames, "last_ns": last_ns, "istrng": istrng,
                "nframe": nframe})
    if decomp:
        dec = parse_decomposition(out)
        if dec:
            (out / "decomp.json").write_text(
                json.dumps(dec, ensure_ascii=False, indent=1), encoding="utf-8")
            res["decomp_file"] = str(out / "decomp.json")
    (out / "mmpbsa.json").write_text(json.dumps(res, ensure_ascii=False, indent=2),
                                     encoding="utf-8")
    if not keep_center and whole.exists():
        whole.unlink()
    return res


def _whole_traj(gmx: str, prod: Path, ndx: Path, out: Path, log: Path) -> Path:
    """把生产轨迹按 -pbc mol -center 处理成整分子、居中，供 MM-PBSA 用。

    原始 md.xtc 里蛋白会跨周期边界，cpptraj 眼里就成了断成两半的分子，
    BOND 项能飙到 1e7 甚至溢出成 '*************'，gmx_MMPBSA 直接报
    “Some energy terms are undefined”。居中后能量项恢复正常量级。
    """
    xtc = prod / "md.xtc"
    dst = out / "center.xtc"
    if dst.exists() and dst.stat().st_mtime >= xtc.stat().st_mtime:
        return dst
    cmd = [gmx, "trjconv", "-s", str(prod / "md.tpr"), "-f", str(xtc),
           "-n", str(ndx), "-o", str(dst), "-pbc", "mol", "-center",
           "-ur", "compact"]
    r = subprocess.run(cmd, input="Protein\nSystem\n", text=True,
                       capture_output=True)
    with open(log, "a", encoding="utf-8") as fh:
        fh.write("\n===== trjconv -pbc mol -center =====\n")
        fh.write((r.stdout or "")[-4000:])
        fh.write((r.stderr or "")[-2000:])
    if r.returncode != 0 or not dst.exists():
        raise RuntimeError(f"trjconv 居中失败，见 {log}")
    return dst


def _n_frames(xtc: Path, gmx: str) -> tuple[int, float]:
    """返回 (帧数, 末帧时间 ps)。"""
    r = subprocess.run([gmx, "check", "-f", str(xtc)], capture_output=True, text=True)
    n, last = 5001, 100000.0
    for line in r.stdout.splitlines():
        s = line.strip()
        if s.startswith("Coords"):
            try:
                n = int(s.split()[1])
            except (IndexError, ValueError):
                pass
        elif s.startswith("Last timestep"):
            try:
                last = float(s.split()[2])
            except (IndexError, ValueError):
                pass
    return n, last


def parse_results(dat: Path) -> dict:
    """从 FINAL_RESULTS.dat 里抓 Δ 项（POISSON BOLTZMANN 段）。"""
    if not dat.exists():
        return {}
    lines = dat.read_text(encoding="utf-8", errors="replace").splitlines()
    out = {}
    started = False
    for line in lines:
        s = line.strip()
        if s.upper().startswith("POISSON BOLTZMANN"):
            started = True
        if not started:
            continue
        if s.startswith("ΔTOTAL") or s.startswith("DELTA TOTAL"):
            q = s.split()
            out["dG"] = float(q[1])
            out["sd"] = float(q[2]) if len(q) > 2 else None
        for t in ("VDWAALS", "EEL", "EPB", "ENPOLAR", "EDISPER", "GGAS", "GSOLV"):
            if s.startswith("Δ" + t):
                q = s.split()
                out[t.lower()] = float(q[1])
    return out


def parse_decomposition(out: Path) -> list[dict]:
    """逐残基分解结果。gmx_MMPBSA 写出 *_DECOMP_MMPBSA.dat。"""
    cands = sorted(out.glob("*DECOMP*.dat")) + sorted(out.glob("*DECOMP*.csv"))
    if not cands:
        return []
    rows = []
    for line in cands[0].read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith(("#", "|", "-", "Residue", "Total")):
            continue
        q = s.replace(",", " ").split()
        if len(q) < 3:
            continue
        try:
            resname = q[1] if not q[1].lstrip("-").isdigit() else "?"
            resid = int(q[1] if q[1].lstrip("-").isdigit() else q[2])
            vals = [float(x) for x in q[-6:]]
        except (ValueError, IndexError):
            continue
        rows.append({"resname": resname, "resid": resid, "dG": vals[-1]})
    return rows


# ---------------------------------------------------------------- 出图
def plot_summary(results: list[dict], outbase: Path, labels=None, want_pdf=True):
    """ΔG_bind 柱状图：每个体系一个柱（多重复取均值），点表示各重复。"""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.family"] = ["Arial", "Microsoft YaHei", "DejaVu Sans"]
    plt.rcParams["axes.linewidth"] = 0.9
    if not results:
        print("没有可画的结果")
        return None
    order = []
    groups: dict[str, list[float]] = {}
    for r in results:
        key = r["job"]
        groups.setdefault(key, []).append(r.get("dG", np.nan))
        if key not in order:
            order.append(key)
    order.sort(key=lambda k: np.nanmean(groups[k]))
    labels = labels or [k.split("__")[-1].replace("_", " ") for k in order]
    means = [np.nanmean(groups[k]) for k in order]
    mins = [np.nanmin(groups[k]) for k in order]
    maxs = [np.nanmax(groups[k]) for k in order]

    fig, ax = plt.subplots(figsize=(max(4.2, 0.85 * len(order) + 2.2), 3.6))
    colors = ["#2b6cb0" if m < 0 else "#c05621" for m in means]
    x = np.arange(len(order))
    ax.bar(x, means, color=colors, width=0.62, alpha=0.85, zorder=2)
    ax.errorbar(x, means, yerr=[np.array(means) - np.array(mins),
                                np.array(maxs) - np.array(means)],
                fmt="none", ecolor="#4a5568", elinewidth=0.9, capsize=3, zorder=3)
    for i, k in enumerate(order):
        vals = [v for v in groups[k] if not np.isnan(v)]
        jitter = np.linspace(-0.09, 0.09, len(vals)) if len(vals) > 1 else [0]
        ax.scatter(np.full(len(vals), x[i]) + jitter, vals, s=18,
                   facecolor="white", edgecolor="#1a202c", linewidth=0.7, zorder=4)
    ax.axhline(0, color="#4a5568", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("$\\Delta G_{bind}$ (kcal mol$^{-1}$)")
    ax.set_title("MM-PBSA binding free energy", fontsize=10.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    plt.tight_layout()
    outbase.parent.mkdir(parents=True, exist_ok=True)
    outs = [outbase.with_suffix(".png")]
    if want_pdf:
        outs.append(outbase.with_suffix(".pdf"))
    for p in outs:
        fig.savefig(p, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  " + ", ".join(p.name for p in outs))
    return outs[0]


def plot_decomposition(rows: list[dict], outbase: Path, title="", top=25,
                       want_pdf=True):
    """逐残基 ΔG 贡献：按 |ΔG| 取前 N 个残基。"""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.family"] = ["Arial", "Microsoft YaHei", "DejaVu Sans"]
    plt.rcParams["axes.linewidth"] = 0.9
    rows = [r for r in rows if abs(r.get("dG", 0)) > 0]
    if not rows:
        print("没有分解结果")
        return None
    rows.sort(key=lambda r: abs(r["dG"]), reverse=True)
    rows = sorted(rows[:top], key=lambda r: r["dG"])
    y = np.arange(len(rows))
    vals = [r["dG"] for r in rows]
    fig, ax = plt.subplots(figsize=(4.6, max(3.0, 0.22 * len(rows) + 1.4)))
    ax.barh(y, vals, color=["#2b6cb0" if v < 0 else "#c05621" for v in vals],
            height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['resname']}{r['resid']}" for r in rows], fontsize=8)
    ax.axvline(0, color="#4a5568", lw=0.8)
    ax.set_xlabel("$\\Delta G$ contribution (kcal mol$^{-1}$)")
    if title:
        ax.set_title(title, fontsize=9.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    plt.tight_layout()
    outbase.parent.mkdir(parents=True, exist_ok=True)
    outs = [outbase.with_suffix(".png")]
    if want_pdf:
        outs.append(outbase.with_suffix(".pdf"))
    for p in outs:
        fig.savefig(p, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  " + ", ".join(p.name for p in outs))
    return outs[0]


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description="MM-PBSA（单作业 / 批量）")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--job", nargs="+", help="流水线作业名")
    g.add_argument("--dir", help="任意作业目录")
    g.add_argument("--plot", nargs="+", metavar="作业", help="只重画已有结果")
    ap.add_argument("--replicates", nargs="+", type=int)
    ap.add_argument("--frames", type=int, default=10, help="抽多少帧（默认 10）")
    ap.add_argument("--last-ns", type=float, default=20.0,
                    help="只用轨迹最后 N ns 取样（默认 20，设 0 表示全程）")
    ap.add_argument("--istrng", type=float, default=0.150)
    ap.add_argument("--keep-center", action="store_true",
                    help="保留居中后的全轨迹 center.xtc（默认算完就删，省空间）")
    ap.add_argument("--decomp", action="store_true",
                    help="逐残基分解（当前不可用，见 README 说明）")
    ap.add_argument("--outdir", help="对比图输出目录，默认 pipeline/figures")
    ap.add_argument("--no-pdf", action="store_true")
    a = ap.parse_args(argv)

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import mdpipe
    cfg = mdpipe.config()

    if a.plot:
        results = []
        for n in a.plot:
            for p in sorted(mdpipe.job_dir(n).glob("08_mmpbsa*/mmpbsa.json")):
                results.append(json.loads(p.read_text(encoding="utf-8")))
        outdir = Path(a.outdir) if a.outdir else mdpipe.PIPE / "figures"
        plot_summary(results, outdir / "MMPBSA_binding", want_pdf=not a.no_pdf)
        return

    if a.decomp:
        print("注意：gmx_MMPBSA 的逐残基分解在这套 GROMACS→ParmEd 拓扑上会因内能项"
              "溢出（'*********'）而失败；本次只算总结合自由能。")
        print("      逐残基相互作用改用 mdpipe.py interactions（见 README）。")

    names = [a.dir] if a.dir else a.job
    results = []
    for n in names:
        job = Path(n) if a.dir else mdpipe.job_dir(n)
        meta = mdpipe.job_meta(job.name) if not a.dir else {}
        nrep = int(meta.get("replicates", 1)) if meta else 1
        reps = a.replicates or list(range(1, nrep + 1))
        for rep in reps:
            print(f"[{job.name}] MM-PBSA 重复 {rep}（{a.frames} 帧）…", flush=True)
            try:
                results.append(run_one(job, rep if nrep > 1 else None,
                                       frames=a.frames, istrng=a.istrng,
                                       decomp=a.decomp,
                                       gmx=cfg.get("gmx") or "gmx",
                                       resname=meta.get("resname", "LIG")
                                       if meta else "LIG",
                                       last_ns=a.last_ns,
                                       keep_center=a.keep_center))
                print(f"  ΔG = {results[-1].get('dG', '?'):.2f} kcal/mol", flush=True)
            except Exception as exc:
                print(f"  失败：{exc}")
                continue
            if a.decomp:
                dp = Path(results[-1].get("decomp_file", ""))
                if dp.exists():
                    rows = json.loads(dp.read_text(encoding="utf-8"))
                    plot_decomposition(rows, Path(results[-1]["dir"]) / "decomp",
                                       title=f"{job.name} 重复 {rep}", want_pdf=not a.no_pdf)

    if results:
        outdir = Path(a.outdir) if a.outdir else mdpipe.PIPE / "figures"
        plot_summary(results, outdir / "MMPBSA_binding", want_pdf=not a.no_pdf)
        csv = outdir / "MMPBSA_binding.csv"
        cols = ["job", "rep", "dG", "sd", "vdwaals", "eel", "epb", "enpolar", "ggas", "gsolv"]
        with open(csv, "w", encoding="utf-8-sig") as fh:
            fh.write(",".join(cols) + "\n")
            for r in results:
                fh.write(",".join(str(r.get(c, "")) for c in cols) + "\n")
        print(f"汇总：{csv}")


if __name__ == "__main__":
    main()
