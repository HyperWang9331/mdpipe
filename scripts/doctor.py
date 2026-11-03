#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""环境自检：装没装、配没配对、能不能出图。

    python scripts/doctor.py
    python scripts/doctor.py --json          # 给安装程序/CI 用

分三档：✓ 可用；! 可选但缺失（功能降级）；✗ 缺了就跑不起来。
只读检测，不改任何文件；每条 ✗ 后面给出下一步怎么做。
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PIPE = ROOT / "pipeline"
CONFIG = PIPE / "config.json"

OK, WARN, MISS = "ok", "warn", "miss"
MARK = {OK: "✓", WARN: "!", MISS: "✗"}

# 安装清单的唯一数据源：文档、doctor --plan、安装说明都照这张表写。
# 级别：必需 = 缺了流程跑不起来；默认路线 = 默认的自动参数化/MM-PBSA 需要，
#       可用自备文件或不做该功能跳过；可选 = 缺了只是功能降级。
# 方式：自动 = 一条 pip/winget 命令，无需交互；半自动 = 要管理员、重启或点安装向导；
#       手动 = 不能也不该由程序代装（驱动、商业/受限软件、账号）。
REQUIREMENTS = [
    ("Python 3.9+", "必需", "半自动",
     "winget install --id Python.Python.3.12 -e（或官网安装包）",
     "安装时勾选 Add python.exe to PATH"),
    ("numpy / scipy / matplotlib", "必需", "自动",
     "pip install numpy scipy matplotlib", "分析、出图的基础"),
    ("GROMACS（CUDA 版）", "必需", "手动",
     "https://manual.gromacs.org", "Windows 无官方 CUDA 二进制；装好把 gmx.exe 填进 config.json"),
    ("NVIDIA 驱动", "建议", "手动",
     "https://www.nvidia.cn/Download/index.aspx",
     "用 GPU 时必须，且版本要满足 GROMACS 的 CUDA 要求；驱动不代装"),
    ("显卡（NVIDIA，建议 ≥6 GB 显存）", "建议", "手动",
     "—", "没有显卡就把 config.json 的 use_gpu 设为 false，纯 CPU 慢约 30 倍"),
    ("WSL2 + Ubuntu", "默认路线", "半自动",
     "管理员 PowerShell：wsl --install -d Ubuntu-24.04",
     "需重启并首次创建用户；配体参数化和 MM-PBSA 需要"),
    ("Miniforge / conda", "默认路线", "半自动",
     "WSL 内下载 Miniforge3 安装脚本并运行", "AmberTools、gmx_MMPBSA 的宿主"),
    ("AmberTools", "默认路线", "半自动",
     "conda create -n amber -c conda-forge ambertools -y",
     "默认配体参数化用；自备 .itp/.gro 可跳过"),
    ("gmx_MMPBSA", "可选", "半自动",
     "conda create -n mmpbsa -c conda-forge ambertools gromacs=2019.6 -y && pip install gmx_MMPBSA",
     "只做 MM-PBSA 时需要"),
    ("pandas", "可选", "自动", "pip install pandas", "表格读写"),
    ("RDKit", "可选", "自动", "pip install rdkit", "2D 相互作用图、糖醛酸判别"),
    ("PyMOL", "可选", "自动", "pip install pymol-open-source", "结构图"),
    ("LigPlot+", "可选", "半自动",
     r"官网下载 zip，解压后把根目录填进 config.json 的 ligplot_home",
     "经典 2D 图；缺失时自动改用 RDKit 图"),
    ("Ghostscript", "可选", "半自动",
     "官网 64 位安装包（或 winget 搜索 ghostscript）", "LigPlot+ 的 .ps 转图"),
    ("OpenBabel", "可选", "半自动",
     "官网安装包或 conda install -c conda-forge openbabel", ".cif 受体、配体加氢"),
    ("Tkinter", "可选", "半自动", "装 Python 时勾选 tcl/tk", "只有图形界面需要"),
]

MODE_TITLES = [
    ("自动", "一、可以自动装（有网络即可，不需要管理员）"),
    ("半自动", "二、半自动（需要管理员、重启或点几下安装向导）"),
    ("手动", "三、只能手动（流程必需，程序不能也不该代装）"),
]


def _plan_text() -> str:
    lines = ["空空如也的电脑：装什么、怎么装",
             "=" * 56,
             "级别：必需 = 缺了跑不起来；默认路线 = 自动参数化/MM-PBSA 需要；"
             "建议 = 影响速度或体验；可选 = 缺了只是功能降级。", ""]
    for mode, title in MODE_TITLES:
        rows = [r for r in REQUIREMENTS if r[2] == mode]
        if not rows:
            continue
        lines.append(title)
        for name, level, _mode, how, note in rows:
            lines.append(f"  · {name}（{level}）")
            lines.append(f"      {how}")
            if note:
                lines.append(f"      说明：{note}")
        lines.append("")
    lines += [
        "【最小可跑组合】",
        "  只分析已有轨迹：Python + numpy/scipy/matplotlib + GROMACS",
        "  从复合物 PDB 跑到 MD：再加 WSL + AmberTools（或自备配体 .itp/.gro）",
        "  加 MM-PBSA：再加 gmx_MMPBSA 环境",
        "  加全套图：再加 RDKit / PyMOL / LigPlot+ / Ghostscript",
        "",
        "“自动”一栏只是 Python 包，前提是机器上已经有 Python 3.9+。",
        "这台机器现在缺什么，直接跑 doctor（不加 --plan）看逐项状态。",
    ]
    return "\n".join(lines)

results: list[dict] = []


def add(name: str, status: str, detail: str, fix: str = "", group: str = ""):
    results.append({"name": name, "status": status, "detail": detail,
                    "fix": fix, "group": group})
    return status


PIP_REQUIRED = ["numpy", "scipy", "matplotlib"]
PIP_OPTIONAL = ["pandas", "rdkit", "pymol-open-source"]
_PIP_MODULE = {"pymol-open-source": "pymol"}


def missing_packages(names) -> list[str]:
    out = []
    for n in names:
        mod = _PIP_MODULE.get(n, n)
        try:
            found = importlib.util.find_spec(mod) is not None
        except (ImportError, ValueError):
            found = False
        if not found:
            out.append(n)
    return out


def pip_command(pkgs) -> list[str]:
    return [sys.executable, "-m", "pip", "install", *pkgs]


def install_pip(names, yes=False) -> int:
    """装缺失的 Python 包。需要网络；不加 --yes 会先列出命令等确认。"""
    pkgs = missing_packages(names)
    if not pkgs:
        print("列出的 Python 包都已安装，不需要动作。")
        return 0
    cmd = pip_command(pkgs)
    print("将执行：" + " ".join(cmd))
    if not yes:
        if not sys.stdin.isatty():
            print("非交互环境，未执行。确认无误后加 --yes 重跑。")
            return 2
        if input("回车开始安装，Ctrl+C 取消：").strip():
            print("已取消。")
            return 2
    return subprocess.run(cmd).returncode


def cfg() -> dict:
    try:
        return json.loads(CONFIG.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _run(cmd, timeout=30, decode_utf16=False):
    try:
        r = subprocess.run(cmd, capture_output=True, timeout=timeout)
    except FileNotFoundError:
        return None
    except subprocess.TimeoutExpired:
        return -1, "", "timeout"
    out = r.stdout
    err = r.stderr
    if decode_utf16 and b"\x00" in out[:40]:
        out = out.decode("utf-16-le", errors="replace")
        err = err.decode("utf-16-le", errors="replace")
    else:
        out = out.decode("utf-8", errors="replace")
        err = err.decode("utf-8", errors="replace")
    return r.returncode, out, err


# ---------------------------------------------------------------- Python
def check_python():
    v = sys.version_info
    version = f"{v.major}.{v.minor}.{v.micro}"
    if v < (3, 9):
        add("Python", MISS, f"{version}（需要 3.9+）",
            "从 https://www.python.org/downloads/ 装 3.11 或更新的版本，"
            "安装时勾选 Add to PATH", "基础")
    else:
        add("Python", OK, f"{version}（{sys.executable}）", group="基础")

    try:
        import tkinter  # noqa: F401
        add("图形界面 tkinter", OK, "可用", group="基础")
    except Exception:
        add("图形界面 tkinter", WARN, "不可用，mdpipe gui 打不开",
            "Windows 版 Python 安装时勾选 tcl/tk；Linux 装 python3-tk", "基础")


PKGS = [
    ("numpy", "必需", "数值计算"),
    ("scipy", "必需", "接触/距离计算"),
    ("matplotlib", "必需", "所有图"),
    ("pandas", "可选", "表格读写"),
    ("rdkit", "可选", "2D 相互作用图、糖醛酸判别"),
]


def check_packages():
    for mod, need, why in PKGS:
        try:
            spec = importlib.util.find_spec(mod)
        except (ImportError, ValueError):
            spec = None
        if spec is not None:
            try:
                from importlib.metadata import version
                ver = version("rdkit" if mod == "rdkit" else mod)
            except Exception:
                ver = "?"
            add(mod, OK, f"{ver}（{why}）", group="Python 包")
        else:
            lvl = MISS if need == "必需" else WARN
            add(mod, lvl, f"未安装（{why}）",
                f"pip install {mod}；本机没有网络时从同版本轮子离线装",
                "Python 包")


# ---------------------------------------------------------------- GROMACS
def check_gmx():
    c = cfg()
    gmx = Path(str(c.get("gmx") or "gmx"))
    src = "config.json"
    if not str(c.get("gmx") or "").strip():
        src = "自动探测"
    if not gmx.exists():
        found = shutil.which(str(gmx))
        if found:
            gmx, src = Path(found), "PATH"
    if not gmx.exists():
        add("GROMACS", MISS, f"找不到 {gmx}（来自 {src}）",
            "运行 mdpipe config --detect --write 自动探测；"
            "也可以把 gmx.exe 的完整路径填进 pipeline/config.json。"
            "Windows 版从 https://manual.gromacs.org 下载（需要 CUDA 版）", "MD 引擎")
        return None
    r = _run([str(gmx), "--version"], timeout=40)
    if not r or r[0] != 0:
        add("GROMACS", MISS, f"{gmx} 无法运行",
            "直接在命令行跑一次 gmx --version 看报错；常见是缺 CUDA 运行库", "MD 引擎")
        return None
    out = r[1]
    ver = cuda = ""
    for line in out.splitlines():
        s = line.strip()
        if s.startswith("GROMACS version:"):
            ver = s.split(":", 1)[1].strip()
        elif s.lower().startswith("gpu support:"):
            cuda = s.split(":", 1)[1].strip()
    detail = f"{ver or '?'} @ {gmx}（{src}）"
    if "CUDA" in cuda.upper():
        add("GROMACS", OK, f"{detail}，GPU 支持：CUDA", group="MD 引擎")
    else:
        add("GROMACS", WARN, f"{detail}，GPU 支持：{cuda or '无'}",
            "当前是 CPU 版 mdrun，200 ns 要跑几周；装 CUDA 版（gmx --version 里应有 CUDA）",
            "MD 引擎")
    return str(gmx)


def check_gpu():
    r = _run(["nvidia-smi",
              "--query-gpu=name,memory.total,driver_version",
              "--format=csv,noheader"])
    if not r or r[0] != 0 or not r[1].strip():
        add("NVIDIA GPU", WARN, "nvidia-smi 不可用",
            "装 NVIDIA 驱动（含 CUDA 支持）；没有独显就把 config.json 的 use_gpu 设为 false"
            "（纯 CPU，慢约 30 倍）", "MD 引擎")
        return
    gpus = [ln.strip() for ln in r[1].splitlines() if ln.strip()]
    add("NVIDIA GPU", OK, "；".join(gpus), group="MD 引擎")


# ---------------------------------------------------------------- WSL / Amber
def check_wsl():
    want = cfg().get("wsl_distro", "Ubuntu-24.04")
    r = _run(["wsl.exe", "-l", "-q"], decode_utf16=True)
    if not r or r[0] != 0:
        add("WSL", WARN, "未安装或无法调用（MM-PBSA / ATB 参数化需要）",
            "管理员 PowerShell 跑 wsl --install -d Ubuntu-24.04", "WSL / Amber")
        return None
    distros = [ln.strip().rstrip("\x00") for ln in r[1].splitlines() if ln.strip()]
    if not any(d.lower() == want.lower() for d in distros):
        add("WSL", WARN, f"没有发行版 {want}（现有：{', '.join(distros) or '无'}）",
            f"wsl --install -d {want}，或把 config.json 的 wsl_distro 改成现有的名字",
            "WSL / Amber")
        return None
    add("WSL", OK, f"{want}", group="WSL / Amber")
    return want


def _wsl_test(distro: str, path: str) -> bool:
    r = _run(["wsl.exe", "-d", distro, "-e", "test", "-x", path])
    return bool(r and r[0] == 0)


def check_amber(distro: str | None):
    c = cfg()
    if not distro:
        add("AmberTools", WARN, "跳过（没有 WSL）", "先装 WSL", "WSL / Amber")
        add("gmx_MMPBSA", WARN, "跳过（没有 WSL）", "先装 WSL", "WSL / Amber")
        return
    amber = c.get("amber_home", "")
    if amber and _wsl_test(distro, amber.rstrip("/") + "/bin/antechamber"):
        add("AmberTools", OK, f"{amber}（config.json: amber_home）", group="WSL / Amber")
    else:
        add("AmberTools", WARN, f"没找到 {amber or '未配置'}",
            "见手册 §14：Miniforge 装到 E 盘后 conda create -n amber -c conda-forge ambertools",
            "WSL / Amber")
    env = c.get("mmpbsa_env") or ""
    if env and _wsl_test(distro, env.rstrip("/") + "/bin/gmx_MMPBSA"):
        add("gmx_MMPBSA", OK, f"{env}（config.json: mmpbsa_env）", group="WSL / Amber")
    else:
        add("gmx_MMPBSA", WARN, f"没找到 {env}",
            "conda create -n mmpbsa -c conda-forge ambertools gromacs=2019.6，"
            "再 pip install gmx_MMPBSA", "WSL / Amber")


# ---------------------------------------------------------------- 出图 / 杂项
def check_figures():
    if importlib.util.find_spec("pymol") is not None:
        add("PyMOL", OK, "可导入（viz 的结构图）", group="出图工具")
    else:
        add("PyMOL", WARN, "未安装，viz 只能出曲线",
            "pip install pymol-open-source；或装官方 PyMOL 并把 scripts 加到 PYTHONPATH",
            "出图工具")
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import interactions
        exedir = interactions.find_ligplot(cfg().get("ligplot_home") or None)
    except Exception:
        exedir = None
    if exedir:
        add("LigPlot+", OK, str(exedir), group="出图工具")
    else:
        add("LigPlot+", WARN, "未检测到（2D 图会退回 RDKit 版本）",
            r"装到 D:\LigPlus\LigPlus，或设置环境变量 LIGPLUS_HOME；"
            r"没有它仍可用 interaction_map", "出图工具")
    gs = shutil.which("gswin64c") or shutil.which("gswin32c") or shutil.which("gs")
    if gs:
        add("Ghostscript", OK, gs, group="出图工具")
    else:
        add("Ghostscript", WARN, "未找到（LigPlot+ 的 .ps 转不了图）",
            "https://www.ghostscript.com/releases/gsdnld.html 装 64 位版", "出图工具")
    ob = shutil.which("obabel")
    if ob:
        add("OpenBabel", OK, ob, group="出图工具")
    else:
        add("OpenBabel", WARN, "未找到（.cif 受体、配体加氢会用到）",
            "conda install -c conda-forge openbabel", "出图工具")


def check_disk():
    try:
        u = shutil.disk_usage(ROOT)
        free = u.free / 2**30
        if free < 10:
            add("磁盘空间", MISS, f"{ROOT.drive} 只剩 {free:.1f} GB",
                "一个 100 ns 体系约 0.6 GB 生产轨迹，先清理或把结果包移到 E 盘",
                "其他")
        elif free < 50:
            add("磁盘空间", WARN, f"{ROOT.drive} 剩 {free:.1f} GB",
                "批量跑之前建议留 50 GB 以上", "其他")
        else:
            add("磁盘空间", OK, f"{ROOT.drive} 剩 {free:.1f} GB", group="其他")
    except Exception:
        pass
    q = PIPE / "queue.json"
    if q.exists():
        add("配置文件", OK, f"{CONFIG.name} / {q.name} 就位", group="其他")
    else:
        add("配置文件", WARN, f"没有 {q}",
            "python scripts/mdpipe.py scan 或 add 之后会生成", "其他")


# ---------------------------------------------------------------- 输出
def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description="MD 流水线环境自检（只读）")
    ap.add_argument("--json", action="store_true", help="输出 JSON")
    ap.add_argument("--report", metavar="PATH",
                    help="同时把报告写进文件（安装包和排查问题时用）")
    ap.add_argument("--plan", action="store_true",
                    help="只打印安装清单（哪些能自动装、哪些必须手动），不检测")
    ap.add_argument("--install-pip", dest="install_pip", action="store_true",
                    help="自动安装缺失的必需 Python 包（numpy/scipy/matplotlib；需网络）")
    ap.add_argument("--with-optional", dest="with_optional", action="store_true",
                    help="配合 --install-pip，一并安装 pandas/rdkit/PyMOL")
    ap.add_argument("--yes", action="store_true", help="不再询问，直接执行")
    ap.add_argument("--no-network", action="store_true", help="跳过 GPU 查询（无显卡机器更快）")
    a = ap.parse_args(argv)

    if a.plan:
        print(_plan_text())
        return 0

    if a.install_pip:
        names = PIP_REQUIRED + (PIP_OPTIONAL if a.with_optional else [])
        rc = install_pip(names, yes=a.yes)
        print("装完重跑一次 doctor 确认。")
        return rc

    check_python()
    check_packages()
    check_gmx()
    if not a.no_network:
        check_gpu()
    distro = check_wsl()
    check_amber(distro)
    check_figures()
    check_disk()

    text = _render(json_out=a.json)
    print(text)
    if a.report:
        path = Path(a.report)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text + "\n", encoding="utf-8")
        print(f"报告已写入 {path}")

    return 1 if any(r["status"] == MISS for r in results) else 0


def _render(json_out: bool = False) -> str:
    if json_out:
        return json.dumps(
            {"results": results,
             "ok": sum(r["status"] == OK for r in results),
             "warn": sum(r["status"] == WARN for r in results),
             "miss": sum(r["status"] == MISS for r in results)},
            ensure_ascii=False, indent=2)
    lines = ["ClawsGO MD 流水线环境检测", "=" * 56]
    group = None
    for r in results:
        if r["group"] != group:
            group = r["group"]
            if group:
                lines.append(f"\n[{group}]")
        lines.append(f" {MARK[r['status']]} {r['name']:<16} {r['detail']}")
        if r["status"] != OK and r["fix"]:
            lines.append(f"   → {r['fix']}")
    n_ok = sum(r["status"] == OK for r in results)
    n_warn = sum(r["status"] == WARN for r in results)
    n_miss = sum(r["status"] == MISS for r in results)
    lines.append("\n" + "-" * 56)
    lines.append(f"可用 {n_ok} 项，可选缺失 {n_warn} 项，必须修复 {n_miss} 项")
    if n_miss:
        lines.append("✗ 的项修好之前，对应命令会直接报错。")
    elif n_warn:
        lines.append("! 的项只会让部分功能降级（例如 2D 图换一种画法），流水线能跑。")
    else:
        lines.append("全部就绪。")
    return "\n".join(lines)


if __name__ == "__main__":
    sys.exit(main())
