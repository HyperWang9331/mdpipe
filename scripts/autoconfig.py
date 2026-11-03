# -*- coding: utf-8 -*-
"""探测本机组件路径，供安装程序和首次运行写入 config.json。

设计原则：
* 只读，不安装、不改系统；
* 探测不到就返回空，由 doctor 逐项提示；
* 所有子进程都有超时，WSL 没装也不会卡住。
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

WIN = sys.platform.startswith("win")

# 可以在 GUI / 命令行里一键重置的路径项
PATH_KEYS = ("gmx", "wsl_home", "amber_home", "mmpbsa_env", "gmx_share",
             "ligplot_home")


def _run(cmd, timeout=12) -> str:
    """跑一个探测命令，返回标准输出；失败返回空串。"""
    try:
        r = subprocess.run(cmd, capture_output=True, timeout=timeout)
    except Exception:
        return ""
    if r.returncode != 0:
        return ""
    raw = r.stdout
    for enc in ("utf-8", "utf-16", "gbk"):
        try:
            text = raw.decode(enc)
        except (UnicodeDecodeError, LookupError):
            continue
        # WSL -l 输出是 UTF-16LE，解码后会带 \x00
        return text.replace("\x00", "")
    return raw.decode("utf-8", errors="replace").replace("\x00", "")


# ---------------------------------------------------------------- GROMACS
def detect_gmx() -> str:
    """按 PATH、常见安装目录、各盘符浅层搜索找 gmx.exe。"""
    found = shutil.which("gmx.exe") or shutil.which("gmx")
    if found:
        return str(Path(found))
    cands: list[Path] = []
    env_home = os.environ.get("GMX_HOME") or os.environ.get("GMXDIR")
    if env_home:
        cands.append(Path(env_home) / "bin" / "gmx.exe")
    for base in (os.environ.get("ProgramFiles", r"C:\Program Files"),
                 os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")):
        cands.extend(sorted(Path(base).glob("GROMACS*/bin/gmx.exe")))
        cands.extend(sorted(Path(base).glob("gromacs*/bin/gmx.exe")))
    for letter in "CDEFGH":
        root = Path(f"{letter}:/")
        if not root.exists():
            continue
        for pat in ("gromacs*/bin/gmx.exe", "GROMACS*/bin/gmx.exe",
                    "gmx*/bin/gmx.exe"):
            try:
                cands.extend(sorted(root.glob(pat)))
            except OSError:
                pass
    cands.append(Path.home() / "gromacs" / "bin" / "gmx.exe")
    for p in cands:
        if p.is_file():
            return str(p)
    return ""


def gmx_share_from(gmx: str) -> str:
    """D:/gromacs2026/bin/gmx.exe -> D:/gromacs2026/share/gromacs/top"""
    if not gmx:
        return ""
    p = Path(gmx)
    top = p.parent.parent / "share" / "gromacs" / "top"
    return top.as_posix() if top.is_dir() else ""


# ---------------------------------------------------------------- WSL
def list_wsl_distros() -> list[str]:
    out = _run(["wsl.exe", "-l", "-q"])
    distros = []
    for line in out.splitlines():
        name = line.strip().strip("*").strip()
        if name and "docker-desktop" not in name.lower():
            distros.append(name)
    return distros


def pick_distro(distros: list[str]) -> str:
    """优先 Ubuntu，其次第一个非 docker 发行版。"""
    for d in distros:
        if d.lower().startswith("ubuntu"):
            return d
    return distros[0] if distros else ""


def wsl_home(distro: str) -> str:
    if not distro:
        return ""
    out = _run(["wsl.exe", "-d", distro, "-e", "sh", "-c", "echo $HOME"])
    home = out.strip().splitlines()[-1].strip() if out.strip() else ""
    return home if home.startswith("/") else ""


def _wsl_exists(distro: str, path: str) -> bool:
    if not distro or not path:
        return False
    r = subprocess.run(["wsl.exe", "-d", distro, "-e", "test", "-x", path],
                       capture_output=True, timeout=15)
    return r.returncode == 0


def _env_candidates(home: str) -> list[str]:
    bases = [home.rstrip("/")] if home else []
    return [f"{b}/{sub}/envs/{name}"
            for b in bases
            for sub in ("miniconda3", "miniforge3", "anaconda3", "mambaforge")
            for name in ("amber", "ambertools", "mmpbsa", "gmxmmpbsa")]


def _wsl_env_dirs(distro: str, exe: str) -> list[str]:
    """扫 WSL 里各盘和家目录下的 conda envs，找带指定可执行文件的环境。"""
    if not distro:
        return []
    script = (
        'for base in /mnt/[c-h]/*/envs /mnt/[c-h]/*/*/envs '
        '"$HOME"/miniconda3/envs "$HOME"/miniforge3/envs "$HOME"/anaconda3/envs; do '
        '[ -d "$base" ] || continue; '
        'for env in "$base"/*; do '
        f'[ -x "$env/bin/{exe}" ] && echo "$env"; '
        'done; done'
    )
    out = _run(["wsl.exe", "-d", distro, "-e", "sh", "-c", script], timeout=25)
    return [ln.strip() for ln in out.splitlines()
            if ln.strip().startswith("/")]


def _dedupe(items):
    seen = set()
    out = []
    for x in items:
        if x and x not in seen:
            seen.add(x)
            out.append(x)
    return out


def detect_amber(distro: str, home: str) -> str:
    cands = ["/opt/amber", "/usr/local/amber", "/opt/conda/envs/amber"]
    cands += [c for c in _env_candidates(home) if "amber" in c.lower()]
    cands += _wsl_env_dirs(distro, "antechamber")
    # 名字完全等于 amber 的优先，其次 ambertools
    cands = sorted(_dedupe(cands),
                   key=lambda p: (p.rstrip("/").rsplit("/", 1)[-1] != "amber", p))
    for p in cands:
        if _wsl_exists(distro, p.rstrip("/") + "/bin/antechamber"):
            return p
    return ""


def detect_mmpbsa(distro: str, home: str) -> str:
    cands = [c for c in _env_candidates(home) if "mmpbsa" in c.lower()]
    cands += _wsl_env_dirs(distro, "gmx_MMPBSA")
    # 名叫 mmpbsa 的环境优先于 gmxmmpbsa 等旧名字
    cands = sorted(_dedupe(cands),
                   key=lambda p: (p.rstrip("/").rsplit("/", 1)[-1] != "mmpbsa", p))
    for p in cands:
        if _wsl_exists(distro, p.rstrip("/") + "/bin/gmx_MMPBSA"):
            return p
    return ""


# ---------------------------------------------------------------- LigPlot+
def detect_ligplot() -> str:
    env = os.environ.get("LIGPLUS_HOME", "").strip()
    cands = [env] if env else []
    cands += [r"D:\LigPlus\LigPlus", r"C:\LigPlus\LigPlus", r"E:\LigPlus\LigPlus"]
    cands.append(str(Path.home() / "LigPlus" / "LigPlus"))
    for c in cands:
        if c and Path(c).is_dir():
            return c
    return ""


# ---------------------------------------------------------------- 汇总
def detect_all(current: dict | None = None) -> dict:
    """返回探测到的路径；探测不到的项不出现。"""
    cur = current or {}
    det: dict[str, str] = {}
    gmx = detect_gmx()
    if gmx:
        det["gmx"] = gmx
        share = gmx_share_from(gmx)
        if share:
            det["gmx_share"] = share
    distros = list_wsl_distros()
    if distros:
        distro = pick_distro(distros)
        det["wsl_distro"] = distro
        home = wsl_home(distro)
        if home:
            det["wsl_home"] = home
        amber = detect_amber(distro, home)
        if amber:
            det["amber_home"] = amber
        mmpbsa = detect_mmpbsa(distro, home)
        if mmpbsa:
            det["mmpbsa_env"] = mmpbsa
    lig = detect_ligplot()
    if lig:
        det["ligplot_home"] = lig
    return det


def merge_detected(cfg: dict, det: dict, overwrite: bool = False) -> dict:
    """把探测结果并进配置；默认只补空项。"""
    out = dict(cfg)
    for key, value in det.items():
        if not value:
            continue
        old = out.get(key, "")
        if overwrite or not old:
            out[key] = value
        elif key in ("gmx", "ligplot_home") and not Path(str(old)).exists():
            out[key] = value
    return out
