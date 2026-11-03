#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_mmpbsa.py —— 一条命令跑完 MM-PBSA（GROMACS → Amber，走 WSL 的 mmpbsa 环境）

为什么要有这个包装：直接调 gmx_MMPBSA 会踩三个坑，脚本把它们按顺序处理掉。
  1. mmpbsa 环境里的 gmx 是 2019 版，读不了 2026 生成的 tpr —— 先用它重建一个结构 tpr；
  2. ParmEd 解析带 #include 的 topol.top 时会把键参数读坏（BOND 变成 5×10^7 kcal/mol）
     —— 先用 2026 的 gmx 做 `grompp -pp` 展开成平铺拓扑再喂给它；
  3. ParmEd 找不到 amber14sb.ff —— 必须给 GMXLIB。

用法：
    python scripts/run_mmpbsa.py --dir MD_FABP5_Glc --frames 6

产物在 MD_FABP5_Glc/08_mmpbsa/ 下，FINAL_RESULTS.dat 里看 ΔTOTAL。
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# 路径优先取环境变量，其次自动探测；都不行就在 main() 里报错并给出设置方法
try:
    import autoconfig
    _DET = autoconfig.detect_all()
except Exception:
    _DET = {}
ENV = os.environ.get("MMPBSA_ENV") or _DET.get("mmpbsa_env", "")
GMXLIB = os.environ.get("GMXLIB") or _DET.get("gmx_share", "")
DISTRO = os.environ.get("WSL_DISTRO") or _DET.get("wsl_distro", "")
WSL_HOME = os.environ.get("WSL_HOME") or _DET.get("wsl_home", "/root")
GMX = os.environ.get("GMX") or _DET.get("gmx", "gmx")


def win_to_wsl(p: Path) -> str:
    s = str(p).replace("\\", "/")
    if len(s) > 1 and s[1] == ":":
        s = "/mnt/" + s[0].lower() + s[2:]
    return s


def run(cmd, cwd=None, log=None, stdin_text=None, check=True, env=None):
    fh = open(log, "w", encoding="utf-8", errors="replace") if log else subprocess.DEVNULL
    try:
        r = subprocess.run(cmd, cwd=str(cwd) if cwd else None, input=stdin_text,
                           stdout=fh, stderr=subprocess.STDOUT, text=True, env=env)
    finally:
        if log:
            fh.close()
    if check and r.returncode != 0:
        raise SystemExit(f"失败: {' '.join(cmd)}\n日志: {log}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="MD 作业目录（含 02_topology/06_prod）")
    ap.add_argument("--gmx", default=GMX, help="gmx.exe 路径（默认自动探测）")
    ap.add_argument("--frames", type=int, default=20, help="取多少帧（均匀抽样）")
    ap.add_argument("--istrng", type=float, default=0.150)
    ap.add_argument("--groups", nargs=2, type=int, default=[1, 13],
                    help="蛋白组与配体组的编号")
    a = ap.parse_args()
    if not ENV or not GMXLIB:
        raise SystemExit(
            "找不到 MM-PBSA 环境。请先设置环境变量：\n"
            "  MMPBSA_ENV=/path/to/envs/mmpbsa\n"
            "  GMXLIB=/path/to/gromacs/share/gromacs/top\n"
            "或在 pipeline/config.json 里填 mmpbsa_env / gmx_share。")

    d = Path(a.dir).resolve()
    top = d / "02_topology"
    prod = d / "06_prod"
    out = d / "08_mmpbsa"
    out.mkdir(exist_ok=True)
    log = out / "run.log"

    # 轨迹帧数：gmx check 末行形如 "Coords          51    20"
    nframe = a.frames
    chk = subprocess.run([a.gmx, "check", "-f", str(prod / "md.xtc")],
                         capture_output=True, text=True).stdout
    for line in chk.splitlines():
        if line.strip().startswith("Coords"):
            try:
                nframe = int(line.split()[1])
            except (IndexError, ValueError):
                pass

    # 1) 展开拓扑（关键：喂给 ParmEd 的必须是平铺文件）
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

    run([a.gmx, "grompp", "-f", str(mdp), "-c", str(top / "system.gro"),
         "-p", str(top / "topol.top"), "-n", str(top / "index.ndx"),
         "-pp", str(out / "processed.top"), "-o", str(out / "pp.tpr"),
         "-maxwarn", "2"], cwd=top, log=out / "grompp_pp.log")

    # 2) 用 mmpbsa 环境的 gmx 重建可读的 tpr
    inner_grompp = (
        f'export PATH={ENV}/bin:$PATH; export GMXLIB={GMXLIB}; '
        f'cd "{win_to_wsl(d)}" && gmx grompp -f 08_mmpbsa/tpr2019.mdp '
        f'-c 02_topology/system.gro -p 02_topology/topol.top -n 02_topology/index.ndx '
        f'-o 08_mmpbsa/complex.tpr -maxwarn 2')
    run(["wsl.exe", "-d", DISTRO, "-e", "env", "-i",
         "PATH=/usr/local/bin:/usr/bin:/bin", "HOME=" + WSL_HOME,
         "bash", "--noprofile", "--norc", "-c", inner_grompp],
        log=out / "grompp2019.log")

    # 3) 参数文件
    (out / "mmpbsa.in").write_text(f"""&general
  startframe=1, endframe={nframe}, interval={max(1, nframe // a.frames)},
  verbose=1, keep_files=0,
/
&pb
  istrng={a.istrng}, fillratio=4.0, radiopt=0,
/
""", encoding="utf-8")

    # 4) 跑 gmx_MMPBSA
    inner = (
        f'export PATH={ENV}/bin:$PATH; export GMXLIB={GMXLIB}; '
        f'cd "{win_to_wsl(d)}/08_mmpbsa" && gmx_MMPBSA -O -i mmpbsa.in '
        f'-cs complex.tpr -ci ../02_topology/index.ndx '
        f'-cg {a.groups[0]} {a.groups[1]} -ct ../06_prod/md.xtc '
        f'-cp processed.top -o FINAL_RESULTS.dat')
    run(["wsl.exe", "-d", DISTRO, "-e", "env", "-i",
         "PATH=/usr/local/bin:/usr/bin:/bin", "HOME=" + WSL_HOME,
         "bash", "--noprofile", "--norc", "-c", inner], log=log)

    res = out / "FINAL_RESULTS.dat"
    if res.exists():
        for line in res.read_text(encoding="utf-8", errors="replace").splitlines():
            if "ΔTOTAL" in line:
                print("MM-PBSA 结合自由能:", line.strip())
                print(f"完整结果: {res}")
                break


if __name__ == "__main__":
    main()
