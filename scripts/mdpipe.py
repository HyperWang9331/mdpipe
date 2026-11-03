#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mdpipe.py —— 蛋白–小分子 MD 批处理流水线（Windows 版）

把《分子动力学操作手册》的十一步串成可断点续跑的作业队列：
投递文件 → 自动建体系 → EM/NVT/NPT/生产 → 汇总指标。
多个作业可以同时排队；准备阶段并行跑，MD 阶段按 GPU 名额排队。

典型用法
--------
    # 1) 投递一个体系（蛋白 PDB + 对接位姿 pdbqt）
    python scripts/mdpipe.py add --protein FABP5.pdb --ligand output.pdbqt --name Gal_1 --charge 0

    # 2) 或者把文件丢进 pipeline/inbox/<作业名>/ 再扫描
    python scripts/mdpipe.py scan
    python scripts/mdpipe.py add --csv jobs.csv      # 一行一个作业

    # 3) 跑
    python scripts/mdpipe.py run --workers 2 --gpu-jobs 2

    # 4) 看状态 / 汇总
    python scripts/mdpipe.py status
    python scripts/mdpipe.py report

目录结构（都在 pipeline/ 下）
--------------------------
    config.json          路径与默认参数（首次运行自动生成）
    queue.json           作业清单与状态
    inbox/<作业名>/       用户投文件的地方
    jobs/<作业名>/        每个作业的 MD_ 目录，结构同手册
    logs/<作业名>/        每一步的完整输出

设计要点
--------
* 每个步骤成功后写 `.done_<step>` 标记，重跑时自动跳过——中断了直接再跑 `run`。
* 准备阶段（拓扑、溶剂化、索引）是 CPU 活，多作业并行；`mdrun` 前必须拿到
  GPU 名额（`--gpu-jobs`），默认 2——本机实测两个 mdrun 同时跑，总吞吐 461 ns/day
  对单个 455–480，基本持平，所以同时跑两个不亏，还能让队列消化得更整齐。
* 所有路径都用绝对路径调用，避免 Windows 上 cd 造成的问题。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = Path(__file__).resolve().parent
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
from split_complex import find_ligand, split_to  # noqa: E402
from selfcheck import check_replicate, job_verdict  # noqa: E402
from notify import notify  # noqa: E402

PIPE = ROOT / "pipeline"
CONFIG_PATH = PIPE / "config.json"
QUEUE_PATH = PIPE / "queue.json"

LOCK = threading.Lock()
# 图形界面点"停止"时置位；每个步骤开始前检查一次，mdrun 之间也检查
STOP = threading.Event()
# --force 时连"日志显示已完成"的保护也绕开，真正从头重跑
FORCE_RERUN = False
# --prep-only：只做到 solvate/ion/index，不碰 GPU（用来验证拓扑，或先准备好再排队）
PREP_ONLY = False
# 同一作业的多个重复是否并行跑。None = 自动：排队作业数少于并行槽位时启用
# （只剩一个体系时，两个重复一起跑能把 GPU 吃满；作业多时并行没有意义）
PARALLEL_REPS = None

DEFAULT_CONFIG = {
    # 组件路径留空，首次运行/安装时由 autoconfig.py 探测；
    # 也可以运行 `mdpipe config --detect --write` 手动写入。
    "gmx": "",
    "wsl_distro": "",
    "amber_home": "",
    # MM-PBSA 用的 conda 环境；gmx_share 留空时从 gmx 路径推出来
    "mmpbsa_env": "",
    "gmx_share": "",
    "wsl_home": "",
    # LigPlot+ 根目录；留空时依次找 D:/LigPlus/LigPlus、C:/LigPlus/... 和 LIGPLUS_HOME
    "ligplot_home": "",
    "gpu_id": 0,
    # false = 纯 CPU 跑 mdrun（没有 NVIDIA 显卡时改这个，慢约 30 倍）
    "use_gpu": True,
    "ntomp": 4,
    "ntomp_em": 8,
    "forcefield": "amber14sb",
    "water": "tip3p",
    "salt_conc": 0.15,
    "box_distance": 1.0,
    "nvt_ps": 500,
    "npt_ps": 1000,
    "prod_ns": 200,
    "gen_seed_base": 1000,
    "gpu_jobs": 2,
    "workers": 2,
    # 轨迹短于这个长度就不算 PCA/自由能景观（短轨迹的 FEL 只是噪声）
    "fel_min_ns": 10,
    # 估算总耗时用：本机两个作业并行时的合计吞吐（ns/day）
    "throughput_ns_per_day": 460,
    # 开工前体检到 FAIL/WARN 时是否自动补 loop / 补侧链
    "auto_repair": True,
    # 队列跑完弹桌面通知 + 蜂鸣；notify_per_job 控制是否每个作业完成都提示
    "notify": True,
    "notify_per_job": False,
}

# 手册里的十个步骤，按顺序执行
PREP_STEPS = ["protein", "pose", "params", "pdb2gmx", "assemble",
              "solvate", "ion", "index"]

# 蛋白力场和配体参数必须同族：GAFF2 配 AMBER，LigParGen 配 OPLS，
# CGenFF 配 CHARMM。混着用的体系能跑完，但参数来自两套不同标定，结果不好解释。
FF_DEFAULT_WATER = {"amber": "tip3p", "opls": "spc", "charmm": "tip3p", "gromos": "spc"}
LIQUID_FF = {"amber": "gaff2", "opls": "opls", "charmm": "cgenff", "gromos": "gromos54a7"}


def ff_family(ff: str) -> str:
    low = ff.lower()
    for fam in ("amber", "opls", "charmm", "gromos"):
        if low.startswith(fam):
            return fam
    return "other"


def gmx_top_dir(cfg) -> Path | None:
    p = Path(cfg.get("gmx", "gmx"))
    try:
        p = p.resolve()
    except OSError:
        return None
    for up in (p.parent.parent, p.parent):
        cand = up / "share" / "gromacs" / "top"
        if cand.is_dir():
            return cand
    return None


def resolve_forcefield(cfg, meta) -> tuple[str, str]:
    """确定本作业用的 (力场, 水模型)，并检查本机确实装了。"""
    ff = (meta.get("forcefield") or cfg.get("forcefield") or "amber14sb").strip()
    water = (meta.get("water") or cfg.get("water") or "").strip()
    explicit_water = bool(meta.get("water"))     # config 里的水模型算"默认值"，可以按力场族换
    fam = ff_family(ff)
    top = gmx_top_dir(cfg)
    if top is not None:
        ffd = top / f"{ff}.ff"
        if not ffd.is_dir():
            avail = sorted(q.name[:-3] for q in top.glob("*.ff"))
            raise RuntimeError(f"本机没有力场 {ff}。可选：{', '.join(avail)}")
        wf = ffd / "watermodels.dat"
        if wf.exists():
            waters = [ln.split()[0] for ln in wf.read_text(errors="replace").splitlines()
                      if ln.strip() and not ln.lstrip().startswith("#")]
            if not water:
                water = FF_DEFAULT_WATER.get(fam, "tip3p")
                if water not in waters:
                    water = waters[0]
            elif water not in waters:
                # config.json 里的水模型是全局默认值。换到别的力场族（比如 gromos54a7
                # 只有 spc/spce）时按族自动配对；只有投递时明确指定了水模型才报错。
                fam_default = FF_DEFAULT_WATER.get(fam, "")
                if not explicit_water and fam_default in waters:
                    water = fam_default
                else:
                    raise RuntimeError(f"{ff} 里没有水模型 {water}。可选：{', '.join(waters)}")
    return ff, water or FF_DEFAULT_WATER.get(fam, "tip3p")


def check_ligand_ff(cfg, meta) -> str:
    """力场和配体参数是否同族。返回警告文本（空串表示没问题）。"""
    ff = meta.get("forcefield") or cfg.get("forcefield") or "amber14sb"
    fam = ff_family(ff)
    ligand = (meta.get("ligand_ff") or "").strip()
    has_itp = bool(meta.get("ligand_itp") or meta.get("ligand_atb"))
    if has_itp:
        return ""
    if fam == "amber":
        return ""
    if ligand and ligand.startswith(LIQUID_FF.get(fam, "?")):
        return ""
    if fam == "gromos":
        return (f"蛋白用 {ff}，配体却走 antechamber/GAFF2 —— 两套参数混用。"
                f"GROMOS 的配体参数要从 ATB 下载（GROMOS 54A7 格式），"
                f"投递时加 --ligand-atb <ATB导出目录>，流水线会自动去重 atomtype。")
    need = LIQUID_FF.get(fam, "对应力场的")
    return (f"蛋白用 {ff}，配体却走 antechamber/GAFF2 —— 两套参数混用。"
            f"要么改回 amber14sb，要么用 {need} 生成配体 itp 后加 "
            f"--ligand-itp/--ligand-gro 投递。")


# --------------------------------------------------------------------------
# 小工具
# --------------------------------------------------------------------------
def load_json(path: Path, default):
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return default


def save_json(path: Path, data):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp.replace(path)


_AUTODETECT_CACHE: dict | None = None


def _autofill_paths(cfg: dict) -> dict:
    """配置里路径为空时探测本机，结果只填进内存，不覆盖用户写的值。"""
    global _AUTODETECT_CACHE
    if cfg.get("gmx") and cfg.get("wsl_distro") and cfg.get("wsl_home"):
        return cfg
    if _AUTODETECT_CACHE is None:
        import autoconfig
        _AUTODETECT_CACHE = autoconfig.detect_all(cfg)
    for key, value in _AUTODETECT_CACHE.items():
        if value and not cfg.get(key):
            cfg[key] = value
    return cfg


def config() -> dict:
    cfg = load_json(CONFIG_PATH, None)
    if cfg is None:
        import autoconfig
        PIPE.mkdir(parents=True, exist_ok=True)
        cfg = autoconfig.merge_detected(dict(DEFAULT_CONFIG), autoconfig.detect_all())
        save_json(CONFIG_PATH, cfg)
    merged = dict(DEFAULT_CONFIG)
    merged.update(cfg)
    return _autofill_paths(merged)


def win_to_wsl(p: Path) -> str:
    """C:\\a\\b -> /mnt/c/a/b"""
    s = str(p).replace("\\", "/")
    if re.match(r"^[A-Za-z]:/", s):
        s = "/mnt/" + s[0].lower() + s[2:]
    return s


def run(cmd, cwd=None, log_path: Path | None = None, stdin_text: str | None = None,
        env: dict | None = None, check=True):
    """执行命令，输出落到日志文件；失败时抛异常并带上日志尾巴。"""
    # gmx 覆盖已有输出时会留 #file.1# 备份，重跑步骤会越堆越多；这里关掉
    if env is None and cmd and str(cmd[0]).lower().endswith(("gmx", "gmx.exe")):
        env = dict(os.environ, GMX_MAXBACKUP="-1")
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        fh = open(log_path, "w", encoding="utf-8", errors="replace")
    else:
        fh = subprocess.DEVNULL
    try:
        r = subprocess.run(cmd, cwd=str(cwd) if cwd else None, input=stdin_text,
                           stdout=fh, stderr=subprocess.STDOUT, text=True,
                           env=env, creationflags=0)
    finally:
        if log_path is not None:
            fh.close()
    if check and r.returncode != 0:
        tail = ""
        if log_path and log_path.exists():
            tail = "\n".join(log_path.read_text(encoding="utf-8", errors="replace")
                             .splitlines()[-12:])
        raise RuntimeError(f"命令失败 ({r.returncode}): {' '.join(map(str, cmd))}\n{tail}")
    return r.returncode


def gmx_cmd(cfg, *args):
    return [cfg["gmx"]] + [str(a) for a in args]


def wsl_sh(cfg, cwd_win: Path, inner: str):
    """在 WSL 的 AmberTools 环境里执行一段 bash。"""
    env_prefix = (f'export AMBERHOME={cfg["amber_home"]}; '
                  f'export PATH=$AMBERHOME/bin:$PATH; ')
    script = f'{env_prefix}cd "{win_to_wsl(cwd_win)}" && {inner}'
    return ["wsl.exe", "-d", cfg["wsl_distro"], "-e", "env", "-i",
            "PATH=/usr/local/bin:/usr/bin:/bin", "HOME=/root",
            "bash", "--noprofile", "--norc", "-c", script]


# --------------------------------------------------------------------------
# 参数溯源：把「这份结果是用什么跑出来的」钉进 job.json
# --------------------------------------------------------------------------
PROV_CACHE: dict = {}


def file_sha256(path: Path, chunk: int = 4 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _gmx_version(cfg) -> str:
    try:
        r = subprocess.run([cfg["gmx"], "--version"], capture_output=True,
                           text=True, timeout=60, errors="replace")
        for line in (r.stdout + r.stderr).splitlines():
            if "GROMACS version" in line:
                return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return "unknown"


def _amber_version(cfg) -> str:
    """AmberTools 是 conda 包，版本号就在 conda-meta 的文件名里。"""
    try:
        cmd = wsl_sh(cfg, ROOT, 'ls "$AMBERHOME/conda-meta/" | grep -i "^ambertools-" | head -1')
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60, errors="replace")
        name = (r.stdout or "").strip().splitlines()
        if name:
            # ambertools-26.0-cuda_None_nompi_py314h...json
            return name[0].split("-")[1]
    except Exception:
        pass
    return "unknown"


def _lib_versions() -> dict:
    out = {}
    for mod in ("numpy", "scipy", "matplotlib", "pandas"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:
            out[mod] = None
    return out


def code_fingerprint() -> dict:
    """脚本内容指纹。改了分析代码，指纹就变，能区分新旧口径的结果。"""
    out = {}
    for fn in ("mdpipe.py", "selfcheck.py", "make_figures.py", "make_fel.py"):
        p = SCRIPTS / fn
        if p.exists():
            out[fn] = file_sha256(p)[:16]
    return out


def provenance(cfg) -> dict:
    """运行环境（力场、gmx、AmberTools、库版本、代码指纹）。缓存，只算一次。"""
    key = (cfg.get("gmx"), cfg.get("amber_home"))
    if key in PROV_CACHE:
        return json.loads(json.dumps(PROV_CACHE[key]))
    ff = cfg.get("forcefield", "")
    ff_hash = None
    top = gmx_top_dir(cfg)
    if top:
        for cand in (top / ff, top / f"{ff}.ff"):
            if (cand / "forcefield.itp").exists():
                ff_hash = file_sha256(cand / "forcefield.itp")[:16]
                break
    import platform
    prov = {
        "recorded": time.strftime("%Y-%m-%d %H:%M:%S"),
        "host": platform.node(),
        "os": platform.platform(),
        "python": platform.python_version(),
        "gmx": _gmx_version(cfg),
        "gmx_path": str(cfg.get("gmx")),
        "amber": _amber_version(cfg),
        "amber_home": str(cfg.get("amber_home")),
        "forcefield": ff,
        "forcefield_itp_sha256": ff_hash,
        "water": cfg.get("water"),
        "libs": _lib_versions(),
        "code": code_fingerprint(),
    }
    PROV_CACHE[key] = prov
    return json.loads(json.dumps(prov))


def job_provenance(cfg, name: str) -> dict:
    """作业级溯源：环境 + 这个作业实际用的 mdp 文件哈希。"""
    prov = provenance(cfg)
    d = job_dir(name)
    mdp = {}
    for p in sorted(d.rglob("*.mdp")):
        mdp[str(p.relative_to(d)).replace("\\", "/")] = file_sha256(p)[:16]
    prov["mdp"] = mdp
    return prov


def cmd_provenance(a):
    names = a.job or (queue()["jobs"] if a.all else [])
    if not names:
        sys.exit("需要 --job 作业名 或 --all")
    cfg = config()
    for n in names:
        meta = job_meta(n)
        if not meta:
            continue
        prov = meta.get("provenance") or {}
        print(f"\n=== {n}")
        if not prov:
            if a.backfill:
                prov = job_provenance(cfg, n)
                prov["backfilled"] = True
                prov["note"] = ("作业跑完后补录；mdp 哈希取自磁盘现存文件，"
                                "代码指纹是补录时的版本")
                meta["provenance"] = prov
                save_meta(n, meta)
                print("  已补录（标记 backfilled）")
            else:
                print("  job.json 里没有溯源信息（作业是在加这个功能之前跑的）")
                print(f"  当前环境：gmx {_gmx_version(cfg)} / AmberTools {_amber_version(cfg)}")
                print("  加 --backfill 可按现存文件补录")
                continue
        for k in ("recorded", "gmx", "amber", "forcefield", "water", "os", "python"):
            if prov.get(k):
                print(f"  {k:<11}{prov[k]}")
        if prov.get("forcefield_itp_sha256"):
            print(f"  {'ff_sha256':<11}{prov['forcefield_itp_sha256']}")
        code = prov.get("code") or {}
        if code:
            print("  code       " + "  ".join(f"{k}:{v}" for k, v in code.items()))
        mdp = prov.get("mdp") or {}
        if mdp:
            print(f"  mdp        {len(mdp)} 个文件")
            for k, v in mdp.items():
                print(f"    {k:<28}{v}")


# 力场族 → mdp 模板子目录。GROMOS 的参数在反应场 + 1.4 nm 截断下拟合，
# 套 amber 那套 PME + 1.0 nm 会让密度和能量系统性偏移，所以单独一套。
MDP_PROFILE = {"gromos": "gromos"}


def tpl_for(cfg, meta, name: str) -> str:
    ff = (meta.get("forcefield") or cfg.get("forcefield") or "").strip()
    sub = MDP_PROFILE.get(ff_family(ff))
    return f"{sub}/{name}" if sub else name


def maxwarn_for(cfg, meta) -> str:
    """grompp 的 -maxwarn 值。

    GROMOS 力场在 Verlet 单截断下会固定报一条"原本是按双截断拟合的"警告，这条是
    力场本身触发的，改 mdp 消除不掉。只对 gromos 族放宽一条，其他力场仍然严格为 0。
    """
    ff = (meta.get("forcefield") or cfg.get("forcefield") or "").strip()
    return "1" if ff_family(ff) == "gromos" else "0"


def mdp_for(cfg, template: str, out: Path, **over):
    """按模板生成 mdp，替换 nsteps 等字段。"""
    src = (ROOT / "mdp_templates" / template).read_text(encoding="utf-8")
    for k, v in over.items():
        key = k.replace("_", "-")
        pattern = re.compile(rf"^{re.escape(key)}\s*=.*$", re.M)
        line = f"{key:<16}= {v}"
        if pattern.search(src):
            src = pattern.sub(line, src, count=1)
        else:
            src += f"\n{line}\n"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(src, encoding="utf-8")


def check_stop():
    if STOP.is_set():
        raise RuntimeError("已停止")


def ps_to_steps(ps: float) -> int:
    return int(round(ps / 0.002))


# --------------------------------------------------------------------------
# 作业模型
# --------------------------------------------------------------------------
def job_dir(name: str, base: Path | None = None) -> Path:
    return (base or PIPE / "jobs") / name


def job_meta(name: str) -> dict:
    return load_json(job_dir(name) / "job.json", {})


def save_meta(name: str, meta: dict):
    save_json(job_dir(name) / "job.json", meta)


def marker(name: str, step: str) -> Path:
    return job_dir(name) / f".done_{step}"


def done(name: str, step: str) -> bool:
    return marker(name, step).exists()


def mark(name: str, step: str):
    marker(name, step).parent.mkdir(parents=True, exist_ok=True)
    marker(name, step).write_text(time.strftime("%Y-%m-%d %H:%M:%S"))


# --------------------------------------------------------------------------
# 步骤实现
# --------------------------------------------------------------------------
def structure_check(path: Path) -> tuple[str, str]:
    """开工前体检，返回 (PASS/WARN/FAIL/?, 一句话说明)。"""
    try:
        import check_protein
        rep = check_protein.analyse(Path(path))
    except Exception as e:                       # 体检本身出错不该挡住作业
        return "?", f"体检没跑成：{e}"
    detail = "；".join(i["what"] for i in rep["issues"][:3]) or "结构完整"
    return rep["verdict"], detail


def repair_structure(cfg, name, src: Path, meta: dict) -> Path:
    """缺 loop / 缺侧链时先补全，再进 pdb2gmx。

    pdb2gmx 遇到主链缺口不会报错：它把链拆成两条、在断口插带电末端，拓扑看着
    正常，体系从第一步就是错的。所以开工前必须先体检。
    """
    d = job_dir(name)
    logs = d / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    verdict, detail = structure_check(src)
    meta["structure_check"] = {"verdict": verdict, "detail": detail,
                               "repaired": False, "source": src.name}
    (logs / "structure_check.log").write_text(
        f"{src}\n结论：{verdict}\n{detail}\n", encoding="utf-8")

    if verdict in ("PASS", "?"):
        return src
    if not cfg.get("auto_repair", True):
        (logs / "structure_check.log").write_text(
            f"{src}\n结论：{verdict}\n{detail}\n\n"
            "未自动补全（config.auto_repair = false）。"
            "跑 `python scripts/prep_receptors_full.py` 补完再投。\n",
            encoding="utf-8")
        return src
    if meta.get("keep_het"):
        (logs / "structure_check.log").write_text(
            f"{src}\n结论：{verdict}\n{detail}\n\n"
            "指定了 keep_het，自动补全会删掉杂原子，这里不补。\n",
            encoding="utf-8")
        return src

    out = d / "00_input" / "_prep"
    cmd = [sys.executable, str(ROOT / "scripts" / "prep_receptors_full.py"),
           str(src), "--out", str(out)]
    for key, flag in (("seqres_dir", "--seqres-dir"), ("seq_dir", "--seq-dir")):
        if meta.get(key):
            cmd += [flag, str(meta[key])]
    run(cmd, log_path=logs / "autoprep.log")

    ready = out / f"{src.stem}_ready.pdb"
    if not ready.exists():
        raise RuntimeError(f"自动补全没产出 {ready}，看 {logs / 'autoprep.log'}")
    meta["structure_check"]["repaired"] = True
    return ready


def step_protein(cfg, name, meta, log):
    d = job_dir(name)
    src = Path(meta["protein"])
    out = d / "00_input" / "protein.pdb"
    (d / "00_input").mkdir(parents=True, exist_ok=True)
    if meta.get("protein_prepared"):
        # 用户自己处理过的文件仍然体检一次，结论记进 job.json，不自动改
        meta["structure_check"] = dict(
            zip(("verdict", "detail"), structure_check(src)), repaired=False,
            source=src.name)
        shutil.copy2(src, out)
        return
    src = repair_structure(cfg, name, src, meta)
    cmd = [sys.executable, str(ROOT / "scripts" / "prep_protein.py"),
           "--in", str(src), "--out", str(out)]
    if meta.get("chain"):
        cmd += ["--chain", meta["chain"]]
    if meta.get("force_altloc"):
        cmd += ["--force-altloc", meta["force_altloc"]]
    if meta.get("keep_het"):
        cmd += ["--keep-het", meta["keep_het"]]
    run(cmd, log_path=log)


def step_pose(cfg, name, meta, log):
    d = job_dir(name)
    (d / "01_ligand").mkdir(parents=True, exist_ok=True)
    pose = Path(meta["ligand"])
    out_pdb = d / "01_ligand" / "pose_mode1.pdb"
    out_mol2 = d / "01_ligand" / "pose_h.mol2"
    if pose.suffix.lower() in (".mol2", ".sdf"):
        # 已经带坐标和氢，直接当配体输入
        shutil.copy2(pose, out_mol2)
        run(["obabel", str(pose), "-O", str(out_pdb)], log_path=log)
        return
    cmd = [sys.executable, str(ROOT / "scripts" / "extract_pose.py"),
           "--pose", str(pose), "--model", str(meta.get("model", 1)),
           "--resname", meta.get("resname", "LIG"),
           "--out", str(out_pdb), "--addh", str(out_mol2)]
    run(cmd, log_path=log)


def _write_posre(itp_path: Path, gro_path: Path, res: str):
    """外部 itp（LigParGen / ATB / CGenFF）通常不带位置约束，按 gro 生成一份。
    assemble_top 把它包在 #ifdef POSRES 里，只有 mdp 明确要求时才生效。"""
    n = 0
    for ln in gro_path.read_text(errors="replace").splitlines()[2:]:
        if len(ln) >= 44:
            n += 1
    lines = [f"; {res} position restraints generated by mdpipe from the gro file",
             "[ position_restraints ]", "; atom  type      fx      fy      fz"]
    lines += [f"{i:6d}     1    1000    1000    1000" for i in range(1, n + 1)]
    itp_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _retag_moleculetype(itp: Path, res: str) -> str | None:
    """外部 itp 里的 moleculetype 名必须和 [molecules] 里的一致，不然 grompp 报错。"""
    text = itp.read_text(errors="replace")
    lines = text.splitlines()
    for i, ln in enumerate(lines):
        if ln.strip().startswith("[") and "moleculetype" in ln:
            for j in range(i + 1, min(i + 6, len(lines))):
                s = lines[j].strip()
                if not s or s.startswith(";"):
                    continue
                name = s.split()[0]
                if name != res:
                    lines[j] = f"{res:<10}" + s[len(name):]
                    itp.write_text("\n".join(lines) + "\n", encoding="utf-8")
                    return name
                return None
    return None


def step_params(cfg, name, meta, log):
    d = job_dir(name)
    res = meta.get("resname", "LIG")
    top = d / "02_topology"
    top.mkdir(parents=True, exist_ok=True)

    # ATB 网页导出的包：先跑 atb_import 整理成流水线能吃的形状（去重 atomtype、
    # 改 moleculetype/残基名、校对原子），再按外部拓扑走
    if meta.get("ligand_atb"):
        src = Path(meta["ligand_atb"])
        if not src.exists():
            raise RuntimeError(f"找不到 ATB 导出包 {src}")
        ff, _ = resolve_forcefield(cfg, meta)
        cmd = [sys.executable, str(ROOT / "scripts" / "atb_import.py"), str(src),
               "--out", str(top), "--resname", res, "--forcefield", ff,
               "--gmx", cfg.get("gmx", "gmx")]
        run(cmd, log_path=top / "atb_import.log")
        meta["ligand_itp"] = str(top / f"{res}.itp")
        meta["ligand_gro"] = str(top / f"{res}.gro")
        meta["ligand_posre"] = str(top / f"posre_{res}.itp")
        with open(log, "a", encoding="utf-8") as fh:
            fh.write(f"配体参数来自 ATB：{src}（力场 {ff}）\n")

    # 外部给的拓扑（LigParGen / ATB / CGenFF），不跑 antechamber
    if meta.get("ligand_itp"):
        itp = Path(meta["ligand_itp"])
        gro = Path(meta.get("ligand_gro") or "")
        if not itp.exists():
            raise RuntimeError(f"找不到配体拓扑 {itp}")
        if not gro or not gro.exists():
            raise RuntimeError("用 --ligand-itp 时必须同时给 --ligand-gro（配体坐标）")
        # ATB 导入器已经直接写进 02_topology，这时源和目标是同一个文件，
        # copy2 在 Windows 上会报 WinError 32（文件被占用）
        def _place(src: Path, dst: Path):
            if src.resolve() != dst.resolve():
                shutil.copy2(src, dst)

        _place(itp, top / f"{res}.itp")
        _place(gro, top / f"{res}.gro")
        old = _retag_moleculetype(top / f"{res}.itp", res)
        log.parent.mkdir(parents=True, exist_ok=True)
        with open(log, "a", encoding="utf-8") as fh:
            fh.write(f"外部配体拓扑：{itp}\n配体坐标：{gro}\n")
            if old:
                fh.write(f"moleculetype {old} 已改名为 {res}\n")
        posre = Path(meta["ligand_posre"]) if meta.get("ligand_posre") else None
        if posre and posre.exists():
            _place(posre, top / f"posre_{res}.itp")
        else:
            _write_posre(top / f"posre_{res}.itp", top / f"{res}.gro", res)
        return

    chg = meta.get("charge", 0)
    inner = (f'bash "{win_to_wsl(ROOT / "scripts" / "param_ligand_wsl.sh")}" '
             f'pose_h.mol2 {res} {chg}')
    run(wsl_sh(cfg, d / "01_ligand", inner), log_path=log)

    acp = d / "01_ligand" / f"{res}.acpype"
    shutil.copy2(acp / f"{res}_GMX.itp", top / f"{res}.itp")
    shutil.copy2(acp / f"posre_{res}.itp", top / f"posre_{res}.itp")
    shutil.copy2(acp / f"{res}_GMX.gro", top / f"{res}.gro")


def step_pdb2gmx(cfg, name, meta, log):
    d = job_dir(name)
    top = d / "02_topology"
    top.mkdir(parents=True, exist_ok=True)
    ff, water = resolve_forcefield(cfg, meta)
    warn = check_ligand_ff(cfg, meta)
    if warn:
        meta.setdefault("warnings", []).append(warn)
    cmd = gmx_cmd(cfg, "pdb2gmx", "-f", d / "00_input" / "protein.pdb",
                  "-o", top / "protein.gro", "-p", top / "topol_protein.top",
                  "-i", top / "posre.itp", "-ff", ff,
                  "-water", water, "-ignh", "-ss")
    run(cmd, cwd=top, log_path=log, stdin_text="y\n")


def step_assemble(cfg, name, meta, log):
    d = job_dir(name)
    top = d / "02_topology"
    res = meta.get("resname", "LIG")
    run([sys.executable, str(ROOT / "scripts" / "assemble_top.py"),
         "--top", str(top / "topol_protein.top"), "--lig-itp", f"{res}.itp",
         "--lig-resname", res, "--lig-count", "1", "--out", str(top / "topol.top")],
        log_path=log)
    run([sys.executable, str(ROOT / "scripts" / "merge_complex.py"),
         "--protein", str(top / "protein.gro"), "--ligand", str(top / f"{res}.gro"),
         "--ligand-resname", res, "--out", str(top / "complex.gro")], log_path=log)


def step_solvate(cfg, name, meta, log):
    d = job_dir(name)
    top = d / "02_topology"
    run(gmx_cmd(cfg, "editconf", "-f", top / "complex.gro", "-o", top / "box.gro",
                "-bt", "dodecahedron", "-d", cfg["box_distance"], "-c"),
        cwd=top, log_path=log)
    run(gmx_cmd(cfg, "solvate", "-cp", top / "box.gro", "-cs", "spc216.gro",
                "-o", top / "solv.gro", "-p", top / "topol.top"), cwd=top, log_path=log)
    mdp = d / "00_mdp" / "em.mdp"
    mdp_for(cfg, tpl_for(cfg, meta, "em.mdp"), mdp)
    run(gmx_cmd(cfg, "grompp", "-f", mdp, "-c", top / "solv.gro", "-p", top / "topol.top",
                "-o", top / "ions.tpr", "-maxwarn", maxwarn_for(cfg, meta)), cwd=top, log_path=log)


def step_ion(cfg, name, meta, log):
    d = job_dir(name)
    top = d / "02_topology"
    run(gmx_cmd(cfg, "genion", "-s", top / "ions.tpr", "-o", top / "system.gro",
                "-p", top / "topol.top", "-pname", "NA", "-nname", "CL",
                "-neutral", "-conc", cfg["salt_conc"]),
        cwd=top, log_path=log, stdin_text="SOL\n")


def step_index(cfg, name, meta, log):
    d = job_dir(name)
    top = d / "02_topology"
    res = meta.get("resname", "LIG")
    run(gmx_cmd(cfg, "make_ndx", "-f", top / "system.gro", "-o", top / "index.ndx"),
        cwd=top, log_path=log, stdin_text=f"1 | r {res}\nq\n")


def _mdrun_prefix(cfg):
    """mdrun 的公共参数。use_gpu=false 时退回纯 CPU（慢约 30 倍，但能跑）。"""
    if not cfg.get("use_gpu", True):
        return [cfg["gmx"], "mdrun", "-ntmpi", "1", "-ntomp", str(cfg["ntomp"])]
    return [cfg["gmx"], "mdrun", "-ntmpi", "1", "-ntomp", str(cfg["ntomp"]),
            "-gpu_id", str(cfg["gpu_id"]), "-nb", "gpu", "-pme", "gpu",
            "-bonded", "gpu"]


def em_step(cfg, name, meta, log, gpu_sem):
    d = job_dir(name)
    top = d / "02_topology"
    em = d / "03_em"
    em.mkdir(exist_ok=True)
    run(gmx_cmd(cfg, "grompp", "-f", d / "00_mdp" / "em.mdp", "-c", top / "system.gro",
                "-p", top / "topol.top", "-n", top / "index.ndx", "-o", em / "em.tpr",
                "-maxwarn", maxwarn_for(cfg, meta)), cwd=em, log_path=log)
    with gpu_sem:
        check_stop()
        cmd = [cfg["gmx"], "mdrun", "-deffnm", "em",
               "-ntmpi", "1", "-ntomp", str(cfg["ntomp_em"])]
        run(cmd, cwd=em, log_path=log)


def _balance(cfg, name, meta, log, gpu_sem, rep: int):
    """NVT → NPT，返回 NPT 目录。"""
    d = job_dir(name)
    top = d / "02_topology"
    sfx = f"_r{rep}" if meta.get("replicates", 1) > 1 else ""
    nvt = d / f"04_nvt{sfx}"
    npt = d / f"05_npt{sfx}"
    nvt.mkdir(exist_ok=True)
    npt.mkdir(exist_ok=True)

    seed = cfg["gen_seed_base"] + rep
    mdp_for(cfg, tpl_for(cfg, meta, "nvt.mdp"), d / "00_mdp" / f"nvt{sfx}.mdp",
            nsteps=ps_to_steps(cfg["nvt_ps"]), gen_seed=seed)
    run(gmx_cmd(cfg, "grompp", "-f", d / "00_mdp" / f"nvt{sfx}.mdp",
                "-c", d / "03_em" / "em.gro", "-r", d / "03_em" / "em.gro",
                "-p", top / "topol.top", "-n", top / "index.ndx",
                "-o", nvt / "nvt.tpr", "-maxwarn", maxwarn_for(cfg, meta)), cwd=nvt, log_path=log)
    with gpu_sem:
        check_stop()
        run(_mdrun_prefix(cfg) + ["-deffnm", "nvt"], cwd=nvt, log_path=log)

    mdp_for(cfg, tpl_for(cfg, meta, "npt.mdp"), d / "00_mdp" / f"npt{sfx}.mdp",
            nsteps=ps_to_steps(cfg["npt_ps"]))
    run(gmx_cmd(cfg, "grompp", "-f", d / "00_mdp" / f"npt{sfx}.mdp",
                "-c", nvt / "nvt.gro", "-r", nvt / "nvt.gro", "-t", nvt / "nvt.cpt",
                "-p", top / "topol.top", "-n", top / "index.ndx",
                "-o", npt / "npt.tpr", "-maxwarn", maxwarn_for(cfg, meta)), cwd=npt, log_path=log)
    with gpu_sem:
        check_stop()
        run(_mdrun_prefix(cfg) + ["-deffnm", "npt"], cwd=npt, log_path=log)
    return npt


def _prod_finished(prod: Path) -> bool:
    """md.log 里出现 'Finished mdrun on' 才算这条轨迹跑完。"""
    log = prod / "md.log"
    if not log.exists():
        return False
    return "Finished mdrun on" in log.read_text(encoding="utf-8", errors="replace")


def prod_step(cfg, name, meta, log, gpu_sem, rep: int):
    d = job_dir(name)
    top = d / "02_topology"
    sfx = f"_r{rep}" if meta.get("replicates", 1) > 1 else ""
    npt = d / f"05_npt{sfx}"
    prod = d / f"06_prod{sfx}"
    prod.mkdir(exist_ok=True)
    if _prod_finished(prod) and not FORCE_RERUN:
        print(f"[{name}] 生产{sfx} 日志显示已完成，跳过重跑", flush=True)
        return prod
    # 2 fs 一步 → 1 ns = 500,000 步
    ns = float(meta.get("prod_ns", cfg["prod_ns"]))
    mdp_for(cfg, tpl_for(cfg, meta, "production.mdp"), d / "00_mdp" / f"production{sfx}.mdp",
            nsteps=int(ns * 500_000))
    run(gmx_cmd(cfg, "grompp", "-f", d / "00_mdp" / f"production{sfx}.mdp",
                "-c", npt / "npt.gro", "-t", npt / "npt.cpt",
                "-p", top / "topol.top", "-n", top / "index.ndx",
                "-o", prod / "md.tpr", "-maxwarn", maxwarn_for(cfg, meta)), cwd=prod, log_path=log)
    with gpu_sem:
        check_stop()
        cmd = _mdrun_prefix(cfg) + ["-deffnm", "md"]
        # 崩溃/断电后重跑：有检查点就从断点接着算，最后一步的不算白跑
        if (prod / "md.cpt").exists():
            print(f"[{name}] 生产{sfx} 从检查点续跑", flush=True)
            cmd += ["-cpi", "md.cpt"]
        run(cmd, cwd=prod, log_path=log)
    return prod


def analysis_step(cfg, name, meta, log, rep: int):
    """轨迹处理 + 六项指标。

    出指标用 md_prot.xtc（Protein_LIG，-pbc mol -center），分子完整、逐帧居中。
    参考结构取它的首帧（ref0.gro），而不是 md.tpr —— 用 tpr 当参考会让配体
    RMSD 出现一个盒子长度的假峰。

    **不要再用 -pbc nojump**：mdrun 写出的 xtc 首帧里就可能有侧链落在相邻镜像
    （Lys40 的 Nζ 离自己的 Cα 6 nm），nojump 是在首帧的基础上做跳变修正，错位
    会被一路带下去。实测这样算出来的骨架 RMSD 从 0.12 nm 假涨到 0.53 nm、
    RMSF 出现 3–6 nm 的尖峰，而 -pbc mol 重拼分子后一切都是正常的。
    """
    d = job_dir(name)
    top = d / "02_topology"
    sfx = f"_r{rep}" if meta.get("replicates", 1) > 1 else ""
    prod = d / f"06_prod{sfx}"
    ana = d / f"07_analysis{sfx}"
    ana.mkdir(exist_ok=True)
    res = meta.get("resname", "LIG")
    tpr, xtc, ndx = prod / "md.tpr", prod / "md.xtc", top / "index.ndx"
    ref, nj = ana / "ref0.gro", ana / "md_prot.xtc"

    def A(tag, args, stdin_text):
        print(f"[{name}] 分析{sfx}: {tag}", flush=True)
        rc = run(gmx_cmd(cfg, *args), cwd=prod, log_path=ana / f"_{tag}.log",
                 stdin_text=stdin_text, check=False)
        if rc != 0:
            print(f"[{name}] 分析{sfx}: {tag} 失败（rc={rc}），见 {ana.name}/_{tag}.log",
                  flush=True)
        return rc

    # 两条轨迹各有用处：
    #   md_center.xtc  全体系 19426 原子，-pbc mol -center → mindist / MM-PBSA
    #   md_prot.xtc    只留 Protein_LIG（2131 原子），-pbc mol -center → 各项指标
    # 指标类工具都要求 -s 的结构与 -f 的轨迹原子数一致；ref0 是从 md_prot 抽的，
    # 所以出指标必须用 md_prot.xtc，拿 md_center.xtc 去配 ref0.gro 会报
    # "Trajectory (19426 atoms) does not match topology (2131 atoms)"。
    A("trjconv_mol",
      ["trjconv", "-s", tpr, "-f", xtc, "-n", ndx,
       "-o", ana / "md_center.xtc", "-pbc", "mol", "-center"], "Protein\nSystem\n")
    A("trjconv_prot",
      ["trjconv", "-s", tpr, "-f", xtc, "-n", ndx,
       "-o", nj, "-pbc", "mol", "-center"], "Protein_LIG\nProtein_LIG\n")
    A("ref0", ["trjconv", "-s", tpr, "-f", nj, "-n", ndx,
               "-o", ref, "-dump", "0"], "Protein_LIG\n")
    if not ref.exists():          # 首帧没出来就退回 tpr 当参考，别把整步卡死
        ref = tpr

    A("rmsd_backbone", ["rms", "-s", ref, "-f", nj, "-n", ndx,
                        "-o", ana / f"rmsd_backbone_{rep}.xvg", "-tu", "ns"],
      "Backbone\nBackbone\n")
    A("rmsd_ligand", ["rms", "-s", ref, "-f", nj, "-n", ndx,
                      "-o", ana / f"rmsd_ligand_{rep}.xvg", "-tu", "ns"],
      f"Protein\n{res}\n")
    A("rmsf", ["rmsf", "-s", ref, "-f", nj, "-n", ndx,
               "-o", ana / f"rmsf_{rep}.xvg", "-res"], "Protein\n")
    A("rg", ["gyrate", "-s", ref, "-f", nj, "-n", ndx,
             "-o", ana / f"rg_{rep}.xvg", "-tu", "ns"], "Protein\n")
    A("sasa", ["sasa", "-s", ref, "-f", nj, "-n", ndx,
               "-o", ana / f"sasa_{rep}.xvg", "-tu", "ns"], "Protein\n")
    # 2026 的 gmx hbond -num 被拆进 hbond-legacy；-s 只认 tpr，-f 用 md_center 轨迹
    A("hbond", ["hbond-legacy", "-s", tpr, "-f", nj, "-n", ndx,
                "-num", ana / f"hbond_{rep}.xvg", "-tu", "ns"],
      f"Protein\n{res}\n")
    A("mindist", ["mindist", "-s", tpr, "-f", ana / "md_center.xtc", "-n", ndx,
                  "-od", ana / "mindist.xvg", "-on", ana / "contacts.xvg"],
      f"{res}\nProtein\n")

    # 自检：用 Rg/SASA/接触数这些与 PBC 无关的量交叉验证 RMSD 类指标
    sc = check_replicate(ana, meta.get("prod_ns"))
    (ana / "selfcheck.json").write_text(
        json.dumps(sc, ensure_ascii=False, indent=2), encoding="utf-8")
    if sc["verdict"] != "ok":
        bad = [c for c in sc["checks"] if c["verdict"] != "ok"]
        print(f"[{name}] 自检{sfx}: {sc['verdict'].upper()} — " + "；".join(
            f"{c['check']}: {c['detail']}" for c in bad), flush=True)
    else:
        print(f"[{name}] 自检{sfx}: 通过", flush=True)
    return ana


ANA_FILES = ["rmsd_backbone_{r}.xvg", "rmsd_ligand_{r}.xvg", "rmsf_{r}.xvg",
             "rg_{r}.xvg", "sasa_{r}.xvg", "hbond_{r}.xvg"]


def analysis_ready(name, meta, rep: int) -> bool:
    """六项指标的 xvg 都在才算分析跑完（老作业缺项时自动补算）。"""
    sfx = f"_r{rep}" if meta.get("replicates", 1) > 1 else ""
    ana = job_dir(name) / f"07_analysis{sfx}"
    return all((ana / f.format(r=rep)).exists() for f in ANA_FILES)


def _traj_ns(ana: Path) -> float:
    """从分析目录的 xvg 读轨迹总时长（ns）。"""
    for pat in ("rmsd_backbone_*.xvg", "rmsd_backbone.xvg"):
        for p in sorted(ana.glob(pat)):
            last = None
            for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
                s = line.strip()
                if s and s[0] not in "#@":
                    last = s
            if last:
                try:
                    return float(last.split()[0])
                except ValueError:
                    pass
    return 0.0


def pca_step(cfg, name, meta, log, rep: int, force: bool = False) -> bool:
    """PCA：协方差矩阵 → 前两个特征向量的二维投影。

    轨迹太短（默认 < 10 ns）时 FEL 只是噪声，跳过；--force 或 fel_min_ns=0 可强制。
    帧数用 -dt 抽稀到千帧量级，200 ns 也不会把协方差矩阵撑爆。
    """
    d = job_dir(name)
    top = d / "02_topology"
    sfx = f"_r{rep}" if meta.get("replicates", 1) > 1 else ""
    ana = d / f"07_analysis{sfx}"
    ref, nj, ndx = ana / "ref0.gro", ana / "md_prot.xtc", top / "index.ndx"
    if not (ref.exists() and nj.exists()):
        print(f"[{name}] PCA 跳过：缺 {ana.name} 的 ref0.gro / md_prot.xtc", flush=True)
        return False
    total = _traj_ns(ana)
    need = float(cfg.get("fel_min_ns", 10))
    if total < need and not force:
        print(f"[{name}] PCA 跳过：轨迹 {total:.2f} ns < {need:g} ns"
              f"（要算就加 --force 或把 config.json 的 fel_min_ns 调小）", flush=True)
        return False

    pca = ana / "pca"
    pca.mkdir(exist_ok=True)
    dt = max(10, int(total * 1000 / 800))     # 目标 ~800 帧
    common = ["-s", ref, "-f", nj, "-n", ndx, "-dt", dt]

    def R(tag, args):
        print(f"[{name}] PCA{sfx}: {tag}", flush=True)
        rc = run(gmx_cmd(cfg, *args), cwd=job_dir(name) / f"06_prod{sfx}",
                 log_path=pca / f"_{tag}.log", stdin_text="Backbone\nBackbone\n",
                 check=False)
        if rc != 0:
            print(f"[{name}] PCA{sfx}: {tag} 失败（rc={rc}）", flush=True)
        return rc

    ok = R("covar", ["covar", *common, "-o", pca / "eigenval.xvg",
                     "-v", pca / "eigenvec.trr", "-av", pca / "avg.pdb",
                     "-l", pca / "covar.log"])
    if ok != 0:
        return False
    R("proj", ["anaeig", *common, "-v", pca / "eigenvec.trr",
               "-eig", pca / "eigenval.xvg", "-proj", pca / "proj.xvg",
               "-first", "1", "-last", "2"])
    R("2d", ["anaeig", *common, "-v", pca / "eigenvec.trr",
             "-eig", pca / "eigenval.xvg", "-2d", pca / "2dproj.xvg",
             "-first", "1", "-last", "2"])
    return (pca / "2dproj.xvg").exists()


def figures_step(cfg, name, meta, log):
    """六联图 + 统计表（+ 有 PCA 数据时的自由能景观）。缺 matplotlib 就跳过，不影响作业。"""
    try:
        sys.path.insert(0, str(ROOT / "scripts"))
        import make_figures
    except Exception as exc:  # noqa: BLE001
        print(f"[{name}] 出图跳过：{exc}", flush=True)
        return None
    print(f"[{name}] 出图 ...", flush=True)
    out = None
    try:
        out = make_figures.plot_job(job_dir(name), name, want_pdf=True)
    except Exception as exc:  # noqa: BLE001
        print(f"[{name}] 出图失败（不影响轨迹与数据）：{exc}", flush=True)
    try:
        import make_fel
        if any((d / "pca" / "2dproj.xvg").exists()
               for d in job_dir(name).glob("07_analysis*")):
            make_fel.plot_job(job_dir(name), name, want_pdf=True)
    except Exception as exc:  # noqa: BLE001
        print(f"[{name}] FEL 出图失败：{exc}", flush=True)
    return out


# --------------------------------------------------------------------------
# 指标
# --------------------------------------------------------------------------
def energy_avg(cfg, edr: Path, term: str) -> float | None:
    """跑 gmx energy，从输出里取该物理量的平均值。"""
    tmp = edr.parent / (edr.stem + f"_{term}.xvg")
    proc = subprocess.run(gmx_cmd(cfg, "energy", "-f", str(edr), "-o", str(tmp)),
                          input=f"{term}\n0\n", text=True,
                          capture_output=True)
    text = proc.stdout
    for line in text.splitlines():
        if line.startswith(term):
            parts = line.split()
            try:
                return float(parts[1])
            except (IndexError, ValueError):
                return None
    return None


def collect_metrics(cfg, name, meta, reps):
    d = job_dir(name)
    m = {}
    em_log = d / "03_em" / "em.log"
    if em_log.exists():
        for line in em_log.read_text(encoding="utf-8", errors="replace").splitlines():
            if "Maximum force" in line:
                m["em_max_force"] = float(line.split("=")[-1].split()[0])
    for rep in reps:
        sfx = f"_r{rep}" if len(reps) > 1 else ""
        nvt = d / f"04_nvt{sfx}"
        npt = d / f"05_npt{sfx}"
        prod = d / f"06_prod{sfx}"
        if (nvt / "nvt.edr").exists():
            m[f"nvt_T{sfx}"] = energy_avg(cfg, nvt / "nvt.edr", "Temperature")
        if (npt / "npt.edr").exists():
            m[f"npt_density{sfx}"] = energy_avg(cfg, npt / "npt.edr", "Density")
        if (prod / "md.log").exists():
            txt = (prod / "md.log").read_text(encoding="utf-8", errors="replace")
            hits = re.findall(r"Performance:\s+([\d.]+)", txt)
            if hits:
                m[f"ns_per_day{sfx}"] = float(hits[-1])
    m["replicates"] = len(reps)
    return m


# --------------------------------------------------------------------------
# 跑一个作业
# --------------------------------------------------------------------------
def run_one_replicate(cfg, name, rep, nreps, gpu_sem, logs):
    """跑一个重复：NVT/NPT → 生产 → 分析 → PCA。

    每个重复的文件都在 *_r{n} 目录里，互不干扰，所以可以被多个线程同时调用；
    共用的准备步骤和 EM 在外面做完再来。
    """
    meta = job_meta(name)
    sfx = f"_r{rep}" if nreps > 1 else ""
    check_stop()
    if not done(name, f"nvt{sfx}"):
        print(f"[{name}] NVT{sfx} ...", flush=True)
        _balance(cfg, name, meta, logs / f"nvt{sfx}.log", gpu_sem, rep)
        mark(name, f"nvt{sfx}")
        mark(name, f"npt{sfx}")
    if not done(name, f"prod{sfx}"):
        print(f"[{name}] 生产{sfx} ...", flush=True)
        prod_step(cfg, name, meta, logs / f"prod{sfx}.log", gpu_sem, rep)
        mark(name, f"prod{sfx}")
    if not done(name, f"ana{sfx}") or not analysis_ready(name, meta, rep):
        print(f"[{name}] 分析{sfx} ...", flush=True)
        analysis_step(cfg, name, meta, logs / f"ana{sfx}.log", rep)
        mark(name, f"ana{sfx}")
    if not done(name, f"pca{sfx}"):
        pca_step(cfg, name, meta, logs / f"pca{sfx}.log", rep)
        mark(name, f"pca{sfx}")


def run_job(cfg, name, gpu_sem, only: list[str] | None = None):
    meta = job_meta(name)
    d = job_dir(name)
    logs = d / "logs"
    logs.mkdir(parents=True, exist_ok=True)

    def logfile(step):
        return logs / f"{step}.log"

    meta["status"] = "running"
    meta["started"] = meta.get("started") or time.strftime("%Y-%m-%d %H:%M:%S")
    meta["provenance"] = provenance(cfg, name)
    save_meta(name, meta)

    steps = [
        ("protein", step_protein), ("pose", step_pose), ("params", step_params),
        ("pdb2gmx", step_pdb2gmx), ("assemble", step_assemble),
        ("solvate", step_solvate), ("ion", step_ion), ("index", step_index),
    ]
    for step, fn in steps:
        check_stop()
        if only and step not in only:
            continue
        if done(name, step):
            continue
        print(f"[{name}] {step} ...", flush=True)
        fn(cfg, name, meta, logfile(step))
        mark(name, step)

    if PREP_ONLY:
        with LOCK:
            meta = job_meta(name)
            meta["status"] = "prepared"
            meta["provenance"] = job_provenance(cfg, name)
            save_meta(name, meta)
        print(f"[{name}] 结构准备完成（--prep-only）：拓扑、溶剂、离子、索引都有了，"
              f"EM 及以后没跑", flush=True)
        return

    reps = list(range(1, meta.get("replicates", 1) + 1))
    # 能量最小化所有重复共用一次
    if not done(name, "em"):
        print(f"[{name}] EM ...", flush=True)
        em_step(cfg, name, meta, logfile("em"), gpu_sem)
        mark(name, "em")

    # 重复默认顺序跑；落单时并行（PARALLEL_REPS），见 cmd_run 的自动判定。
    # 准备步骤和 EM 在循环外只做一次，并行部分只碰各自的 *_r{n} 目录。
    npar = 1
    if len(reps) > 1 and PARALLEL_REPS:
        npar = min(len(reps), max(1, int(cfg.get("gpu_jobs", 2))))
    if npar > 1:
        print(f"[{name}] {len(reps)} 个重复并行跑（同时 {npar} 个）", flush=True)
        with ThreadPoolExecutor(max_workers=npar) as pool:
            futs = [pool.submit(run_one_replicate, cfg, name, rep, len(reps), gpu_sem, logs)
                    for rep in reps]
            for f in as_completed(futs):
                f.result()          # 某个重复失败时照常抛给 worker
    else:
        for rep in reps:
            run_one_replicate(cfg, name, rep, len(reps), gpu_sem, logs)

    if not done(name, "fig"):
        figures_step(cfg, name, meta, logfile("fig"))
        mark(name, "fig")

    with LOCK:
        meta = job_meta(name)
        meta["status"] = "done"
        meta["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
        meta["metrics"] = collect_metrics(cfg, name, meta, reps)
        # 这时 mdp 都已写完，补上本作业实际用到的参数文件指纹与代码指纹
        meta["provenance"] = job_provenance(cfg, name)
        save_meta(name, meta)
    print(f"[{name}] 完成", flush=True)


def worker(cfg, name, gpu_sem, only=None):
    try:
        run_job(cfg, name, gpu_sem, only)
        if cfg.get("notify_per_job"):
            notify(f"作业完成：{short_label(name)}",
                   f"{name} 已跑完 {job_meta(name).get('prod_ns', '?')} ns", "info",
                   sound=False)
        return name, True, ""
    except Exception as exc:  # noqa: BLE001
        stopped = STOP.is_set()
        with LOCK:
            meta = job_meta(name)
            meta["status"] = "stopped" if stopped else "failed"
            meta["error"] = str(exc)[-2000:]
            save_meta(name, meta)
        print(f"[{name}] {'已停止' if stopped else '失败'}: {exc}", flush=True)
        if cfg.get("notify_per_job") and not stopped:
            notify(f"作业失败：{short_label(name)}", str(exc)[-200:], "fail")
        return name, False, str(exc)


# --------------------------------------------------------------------------
# 队列维护
# --------------------------------------------------------------------------
def queue() -> dict:
    return load_json(QUEUE_PATH, {"jobs": []})


def job_steps(meta: dict) -> list[str]:
    """一个作业会留下的全部步骤标记，用来算进度。"""
    reps = int(meta.get("replicates", 1))
    steps = list(PREP_STEPS) + ["em"]
    for r in range(1, reps + 1):
        sfx = f"_r{r}" if reps > 1 else ""
        steps += [f"nvt{sfx}", f"npt{sfx}", f"prod{sfx}", f"ana{sfx}", f"pca{sfx}"]
    steps.append("fig")
    return steps


def progress(names=None) -> dict:
    """整体进度与剩余时间估算。

    剩余 ns 只算还没跑生产的那部分；NVT/NPT 的固定开销按每条重复 ~1.7 ns 折算。
    吞吐用 config 的 throughput_ns_per_day（默认 460，是本机两个作业并行的合计实测值），
    有实测 ns/day 时取中位数。
    """
    cfg = config()
    q = queue()
    names = names or q["jobs"]
    total = done_n = jobs_done = 0
    remain_ns = 0.0
    speeds = []
    for n in names:
        meta = job_meta(n)
        if not meta:
            continue
        steps = job_steps(meta)
        total += len(steps)
        done_n += sum(1 for s in steps if done(n, s))
        reps = int(meta.get("replicates", 1))
        finished = meta.get("status") == "done"
        if finished:
            jobs_done += 1
        spd = (meta.get("metrics") or {}).get("ns_per_day")
        if isinstance(spd, (int, float)) and spd > 0:
            speeds.append(spd)
        if not finished:
            for r in range(1, reps + 1):
                sfx = f"_r{r}" if reps > 1 else ""
                if not done(n, f"prod{sfx}"):
                    remain_ns += float(meta.get("prod_ns", cfg["prod_ns"]))
                if not done(n, f"npt{sfx}"):
                    remain_ns += 1.7          # NVT 0.5 + NPT 1 + EM 的折算
    speed = statistics.median(speeds) if speeds else 0.0
    if speed <= 0:
        speed = float(cfg.get("throughput_ns_per_day", 460))
    eta_h = remain_ns / speed * 24 if remain_ns > 0 else 0.0
    return {
        "total_steps": total, "done_steps": done_n,
        "pct": (100.0 * done_n / total) if total else 0.0,
        "jobs_total": len(names), "jobs_done": jobs_done,
        "remaining_ns": remain_ns, "eta_hours": eta_h, "speed": speed,
    }


def save_queue(q):
    PIPE.mkdir(parents=True, exist_ok=True)
    save_json(QUEUE_PATH, q)


def add_job(meta: dict) -> bool:
    name = meta["name"]
    d = job_dir(name)
    if d.exists() and (d / "job.json").exists():
        print(f"作业 {name} 已存在，跳过（要重来请先 mdpipe.py remove {name}）")
        return False
    for sub in ("00_input", "00_mdp", "01_ligand", "02_topology",
                "03_em", "04_nvt", "05_npt", "06_prod", "07_analysis", "logs"):
        (d / sub).mkdir(parents=True, exist_ok=True)
    meta.setdefault("resname", "LIG")
    meta.setdefault("model", 1)
    meta.setdefault("replicates", 1)
    meta.setdefault("status", "pending")
    meta["created"] = time.strftime("%Y-%m-%d %H:%M:%S")
    save_meta(name, meta)
    q = queue()
    if name not in q["jobs"]:
        q["jobs"].append(name)
        save_queue(q)
    print(f"已登记 {name}: {Path(meta['protein']).name} + {Path(meta['ligand']).name}")
    return True


def detect_charge(cfg, protein, ligand) -> int | None:
    """读配体的形式电荷之和。

    做法：pdbqt 先用 extract_pose 取第一个构象并补氢成 mol2；再用
    `obabel x.mol2 -oreport` 的 "TOTAL CHARGE" 行。只要装了 OpenBabel 就能用，
    不依赖 RDKit（conda 基础环境里没有 RDKit，也能正确判定）。
    """
    lig = Path(ligand)
    work = PIPE / "_charge_probe"
    work.mkdir(parents=True, exist_ok=True)
    try:
        if lig.suffix.lower() in (".mol2", ".sdf"):
            target = lig
        else:
            target = work / "pose.mol2"
            run([sys.executable, str(ROOT / "scripts" / "extract_pose.py"),
                 "--pose", str(lig), "--model", "1",
                 "--out", str(work / "pose.pdb"), "--addh", str(target)],
                log_path=work / "probe.log")
        r = subprocess.run(["obabel", str(target), "-oreport"],
                           capture_output=True, text=True)
        out = r.stdout or ""
        if r.returncode != 0 or "FILENAME" not in out:
            return None
        # 中性分子 obabel 不打印 TOTAL CHARGE 行，缺省即 0
        m = re.search(r"TOTAL CHARGE:\s*(-?\d+)", out)
        return int(m.group(1)) if m else 0
    except Exception:
        return None
    finally:
        shutil.rmtree(work, ignore_errors=True)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
LIGANT_EXTS = (".pdbqt", ".mol2", ".sdf")
PROTEIN_EXTS = (".pdb", ".ent")


def sanitize_name(s: str) -> str:
    """作业名只用字母数字、下划线、连字符和中文，去掉空格与路径符号。"""
    s = re.sub(r"[^\w一-鿿.-]+", "_", str(s)).strip("_.")
    return s or "job"


def find_files(paths) -> list[Path]:
    """把文件/目录混合列表展开成具体文件，按扩展名筛选。"""
    out: list[Path] = []
    for p in paths or []:
        p = Path(p)
        if p.is_dir():
            for f in sorted(p.iterdir()):
                if f.is_file() and f.suffix.lower() in PROTEIN_EXTS + LIGANT_EXTS:
                    out.append(f)
        elif p.is_file():
            out.append(p)
        else:
            print(f"找不到 {p}")
    return out


def detect_complex(path) -> str | None:
    """这个 PDB 看着像"蛋白 + 配体"复合物时，返回配体残基名；否则 None。"""
    p = Path(path)
    if p.suffix.lower() not in PROTEIN_EXTS:
        return None
    try:
        lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
        lig, want_atom = find_ligand(lines)
    except Exception:
        return None
    n = sum(1 for l in lines if l[17:20].strip().upper() == lig
            and (l.startswith("HETATM") or (want_atom and l.startswith("ATOM"))))
    return lig if n >= 6 else None      # 一两个原子多半是离子、不是配体


def split_complex_job(src, jobname: str, resname: str | None = None):
    """把复合物拆到 pipeline/splits/<作业名>/，返回 (protein.pdb, ligand.mol2)。"""
    d = PIPE / "splits" / sanitize_name(jobname)
    return split_to(Path(src), d, resname)


def plan_pairs(pairs, base_name: str | None = None) -> list[tuple[str, Path, Path]]:
    """把"蛋白–配体对"清单变成 (作业名, 蛋白, 配体)。

    - 只有一对：作业名 = base_name 或配体文件名
    - 同一蛋白配多个配体：作业名 = 配体文件名（有 base_name 时加前缀）
    - 多个蛋白（任意组合）：作业名 = 蛋白名__配体名
    重名自动加 _2、_3。
    """
    n = len(pairs)
    same_protein = len({str(p.resolve()) for p, _ in pairs}) == 1
    plan, seen = [], {}
    for p, l in pairs:
        if base_name and n == 1:
            raw = base_name
        elif base_name:
            raw = f"{p.stem}__{base_name}__{l.stem}" if not same_protein else f"{base_name}__{l.stem}"
        elif n == 1:
            raw = l.stem
        elif same_protein:
            raw = l.stem
        else:
            raw = f"{p.stem}__{l.stem}"
        name = sanitize_name(raw)
        if name in seen:
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 1
        plan.append((name, p, l))
    return plan


def plan_jobs(proteins, ligands, base_name: str | None = None) -> list[tuple[str, Path, Path]]:
    """全组合（笛卡尔积）：每个蛋白配每个配体。"""
    return plan_pairs([(p, l) for p in proteins for l in ligands], base_name)


def resolve_token(tok: str, exts, search_dirs) -> Path | None:
    """把一个名字解析成文件：先按路径找，再到给定目录里按文件名或去掉扩展名找。"""
    t = tok.strip().strip('"').strip("'")
    if not t:
        return None
    p = Path(t)
    if p.is_file():
        return p
    low = t.lower()
    for d in search_dirs or []:
        if not d:
            continue
        d = Path(d)
        if not d.is_dir():
            continue
        for f in sorted(d.iterdir()):
            if f.is_file() and f.suffix.lower() in exts and \
                    (f.name.lower() == low or f.stem.lower() == low):
                return f
    return None


def parse_pair_specs(specs, protein_dirs=(), ligand_dirs=()) -> list[tuple[Path, Path]]:
    """解析 `蛋白::配体1,配体2` 形式的配对（用 `::` 分隔，避免和盘符的冒号打架）。

    例：`FABP5.pdb::Glc,Gal` → (FABP5, Glc) + (FABP5, Gal)
    蛋白和配体都可以给全路径。
    """
    protein_dirs = list(protein_dirs or [])
    ligand_dirs = list(ligand_dirs or [])
    pairs, bad, seen = [], [], set()
    for spec in specs:
        if "::" not in spec:
            bad.append(spec)
            continue
        ptok, rest = spec.split("::", 1)
        prot = resolve_token(ptok, PROTEIN_EXTS, protein_dirs)
        if prot is None:
            bad.append(f"{ptok}（蛋白）")
            continue
        for lt in rest.split(","):
            lt = lt.strip()
            if not lt:
                continue
            lig = resolve_token(lt, LIGANT_EXTS, ligand_dirs + [prot.parent])
            if lig is None:
                bad.append(f"{lt}（配体）")
                continue
            key = (str(prot.resolve()), str(lig.resolve()))
            if key not in seen:
                seen.add(key)
                pairs.append((prot, lig))
    for b in bad:
        print(f"  解析不了：{b}")
    return pairs


def load_pairs_table(path: Path) -> list[dict]:
    """读配对表：CSV/TSV，必须有 protein 和 ligand 两列。

    可选列：name / charge / replicates / ns / resname / model。
    """
    import csv
    text = Path(path).read_text(encoding="utf-8-sig")
    first = text.splitlines()[0] if text else ""
    delim = "\t" if "\t" in first else (";" if ";" in first and "," not in first else ",")
    rows = []
    for row in csv.DictReader(text.splitlines(), delimiter=delim):
        row = {(k or "").strip().lower(): (v or "").strip()
               for k, v in row.items() if k}
        if row.get("protein") and row.get("ligand"):
            rows.append(row)
    return rows


def resolve_table_rows(rows, base: Path) -> list[dict]:
    """配对表里的相对路径按表所在目录解析成绝对路径。"""
    out = []
    for row in rows:
        row = dict(row)
        for k in ("protein", "ligand"):
            p = Path(row[k])
            row[k] = str((p if p.is_absolute() else base / p).resolve())
        out.append(row)
    return out


def _register_plan(cfg, plan, a, extras=None):
    """把 (作业名, 蛋白, 配体) 清单逐个登记。extras 按作业名附加字段。"""
    extras = extras or {}
    if len(plan) > 50 and not a.yes:
        print(f"这次会登记 {len(plan)} 个作业。确认无误请加 --yes。")
        return 0
    print(f"将登记 {len(plan)} 个作业")
    added = 0
    for name, prot, lig in plan:
        charge = a.charge
        if charge is None:
            charge = detect_charge(cfg, prot, lig)
            if charge is None:
                print(f"  {name}: 未能自动判定净电荷，按 0 处理——请核对 job.json")
                charge = 0
            else:
                print(f"  {name}: 净电荷 {charge}")
        meta = {
            "name": name, "protein": str(Path(prot).resolve()),
            "ligand": str(Path(lig).resolve()), "charge": int(charge),
            "resname": a.resname, "model": a.model, "replicates": a.replicates,
            "prod_ns": a.ns, "chain": a.chain, "force_altloc": a.force_altloc,
            "keep_het": a.keep_het, "protein_prepared": a.prepared,
            "seq_dir": a.seq_dir, "seqres_dir": a.seqres_dir,
            "forcefield": a.forcefield, "water": a.water,
            "ligand_itp": a.ligand_itp, "ligand_gro": a.ligand_gro,
            "ligand_posre": a.ligand_posre, "ligand_ff": a.ligand_ff,
            "ligand_atb": a.ligand_atb,
        }
        meta.update(extras.get(name, {}))
        warn = check_ligand_ff(cfg, meta)
        if warn:
            print(f"  {name}: 注意 —— {warn}")
        if add_job(meta):
            added += 1
    return added


def _plan_from_complex(a):
    """把复合物 PDB 拆开，生成作业清单。返回 (plan, extras)。"""
    plan, extras = [], {}
    files = [Path(p) for p in a.complex]
    seen = {}
    for src in files:
        name = sanitize_name(a.name or src.stem)
        if a.name and len(files) > 1:
            name = sanitize_name(f"{a.name}__{src.stem}")
        if name in seen:                    # 两个复合物同名（不同目录）
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 1
        try:
            prot, lig = split_complex_job(src, name, a.lig_resname)
        except Exception as e:
            print(f"  {src.name}: 拆分失败——{e}")
            continue
        print(f"  {src.name}: 已拆出 {prot.name} + {lig.name}")
        plan.append((name, prot, lig))
        extras[name] = {"source_complex": str(src.resolve()),
                        "complex_resname": a.lig_resname}
    return plan, extras


def cmd_complex(a):
    """复合物专用入口：只吃「蛋白 + 配体」的复合物 PDB，自动拆分登记。"""
    if a.ns is None:
        a.ns = config()["prod_ns"]
    a.complex = a.files
    a.protein = a.protein_dir = a.ligand = a.ligand_dir = None
    a.pair = a.csv = None
    cmd_add(a)


def cmd_add(a):
    cfg = config()
    PIPE.mkdir(parents=True, exist_ok=True)
    if a.csv:
        table = Path(a.csv).resolve()
        for row in resolve_table_rows(load_pairs_table(table), table.parent):
            for k in ("protein", "ligand"):
                if not Path(row[k]).is_file():
                    print(f"  找不到 {row[k]}")
            meta = dict(row)
            meta["charge"] = int(row.get("charge") or 0)
            meta["replicates"] = int(row.get("replicates") or 1)
            meta["prod_ns"] = float(row.get("ns") or row.get("prod_ns") or cfg["prod_ns"])
            meta["name"] = meta.get("name") or sanitize_name(
                f"{Path(row['protein']).stem}__{Path(row['ligand']).stem}")
            add_job(meta)
        return

    # 0) 直接给复合物：自动拆
    if a.complex:
        plan, extras = _plan_from_complex(a)
        if plan:
            _register_plan(cfg, plan, a, extras)
        return

    protein_dirs = [Path(a.protein_dir)] if a.protein_dir else []
    ligand_dirs = [Path(a.ligand_dir)] if a.ligand_dir else []

    # 1) 显式配对：--pair "P.pdb::a,b" （可重复）
    if a.pair:
        pairs = parse_pair_specs(a.pair, protein_dirs, ligand_dirs)
        if not pairs:
            print("没有解析出任何配对；格式是 --pair \"蛋白::配体1,配体2\"")
            return
        plan = plan_pairs(pairs, a.name)
    else:
        proteins = find_files([a.protein] if a.protein else [])
        if not proteins and a.protein_dir:
            proteins = find_files([a.protein_dir])
        proteins = [p for p in proteins if p.suffix.lower() in PROTEIN_EXTS]
        if a.ligand:
            ligands = find_files(a.ligand)
        elif a.ligand_dir:
            ligands = find_files([a.ligand_dir])
        else:
            ligands = []
        ligands = [l for l in ligands if l.suffix.lower() in LIGANT_EXTS]

        # 只给了蛋白：如果它其实是复合物，自动拆开继续
        if proteins and not ligands and len(proteins) == 1:
            lig_res = detect_complex(proteins[0])
            if lig_res:
                name = sanitize_name(a.name or proteins[0].stem)
                try:
                    prot, ligf = split_complex_job(proteins[0], name, None)
                except Exception as e:
                    print(f"拆分失败——{e}")
                    return
                print(f"{proteins[0].name} 里有配体残基 {lig_res}，按复合物拆开")
                _register_plan(cfg, [(name, prot, ligf)], a,
                               {name: {"source_complex": str(proteins[0].resolve()),
                                       "complex_resname": lig_res}})
                return

        if not proteins or not ligands:
            print("缺少蛋白或配体文件；用 --protein/--protein-dir 与 --ligand/--ligand-dir 指定，"
                  "或用 --complex（复合物 PDB）、--pair、--csv")
            return
        plan = plan_jobs(proteins, ligands, a.name)

    _register_plan(cfg, plan, a)


def cmd_scan(a):
    cfg = config()
    inbox = PIPE / "inbox"
    if not inbox.exists():
        print(f"{inbox} 不存在；先建目录并放文件")
        return
    for sub in sorted(p for p in inbox.iterdir() if p.is_dir()):
        pdbs = [f for f in sub.iterdir() if f.suffix.lower() in PROTEIN_EXTS]
        ligs = [f for f in sub.iterdir() if f.suffix.lower() in LIGANT_EXTS]
        # 只放了一个复合物 PDB：自动拆开，当成一个完整体系
        if pdbs and not ligs:
            lig_res = detect_complex(pdbs[0]) if len(pdbs) == 1 else None
            if lig_res:
                try:
                    prot, ligf = split_complex_job(pdbs[0], sub.name, None)
                except Exception as e:
                    print(f"跳过 {sub.name}：拆分失败——{e}")
                    continue
                print(f"{sub.name}: {pdbs[0].name} 是复合物（配体 {lig_res}），已拆开")
                charge = None
                cfile = sub / "charge.txt"
                if cfile.exists():
                    charge = int(cfile.read_text().strip())
                if charge is None:
                    charge = detect_charge(cfg, prot, ligf)
                add_job({"name": sub.name, "protein": str(prot.resolve()),
                         "ligand": str(ligf.resolve()), "charge": int(charge or 0),
                         "resname": "LIG", "model": 1, "replicates": 1,
                         "prod_ns": cfg["prod_ns"],
                         "source_complex": str(pdbs[0].resolve())})
                continue
            print(f"跳过 {sub.name}：只有一个 PDB 又没有可识别的配体；"
                  f"如果它其实是复合物，检查里面配体的残基名")
            continue
        if not pdbs or not ligs:
            print(f"跳过 {sub.name}：需要至少 1 个蛋白 (.pdb) 和 1 个配体 (.pdbqt/.mol2/.sdf)")
            continue
        # 一对多：一个蛋白 + 多个配体，每个配体一个作业；多对多做交叉组合
        base = sub.name if (len(pdbs) > 1 or len(ligs) > 1) else None
        plan = plan_jobs(sorted(pdbs), sorted(ligs), base)
        cfile = sub / "charge.txt"
        fixed_charge = int(cfile.read_text().strip()) if cfile.exists() else None
        for name, prot, lig in plan:
            charge = fixed_charge
            if charge is None:
                charge = detect_charge(cfg, prot, lig)
            meta = {"name": name, "protein": str(prot.resolve()),
                    "ligand": str(lig.resolve()), "charge": int(charge or 0),
                    "resname": "LIG", "model": 1, "replicates": 1,
                    "prod_ns": cfg["prod_ns"]}
            add_job(meta)


def cmd_run(a):
    cfg = config()
    global PREP_ONLY
    PREP_ONLY = bool(getattr(a, "prep_only", False))
    for k, v in (("nvt_ps", a.nvt_ps), ("npt_ps", a.npt_ps), ("prod_ns", a.ns)):
        if v:
            cfg[k] = v
    q = queue()
    names = a.job or q["jobs"]
    names = [n for n in names if job_meta(n)]
    if a.only_failed:
        names = [n for n in names if job_meta(n).get("status") == "failed"]
    # 已完成但还没出图的作业也排进来：所有步骤都有标记，只会补跑最后的出图
    todo = [n for n in names
            if job_meta(n).get("status") not in ("done",) or a.force
            or not done(n, "fig")]
    # 生产时长存在每个作业的 job.json 里，--ns 要写回去才生效
    if a.ns:
        for n in todo:
            m = job_meta(n)
            m["prod_ns"] = float(a.ns)
            save_meta(n, m)
    if a.force:
        global FORCE_RERUN
        FORCE_RERUN = True
        for n in todo:
            for mk in job_dir(n).glob(".done_*"):
                mk.unlink()
    if not todo:
        print("没有待跑作业")
        return
    workers = a.workers or cfg.get("workers", 2)
    gpu_jobs = a.gpu_jobs or cfg.get("gpu_jobs", 2)
    # 落单：待跑作业比并行槽位少时，把每个作业的重复并行跑，把 GPU 吃满。
    # 显式 --parallel-replicates / --serial-replicates 优先于自动判定。
    global PARALLEL_REPS
    explicit = getattr(a, "parallel_replicates", None)
    PARALLEL_REPS = bool(explicit) if explicit is not None else (len(todo) < workers)
    print(f"共 {len(todo)} 个作业，并行 {workers}，GPU 同时跑 {gpu_jobs} 个"
          + ("；落单作业的重复将并行跑" if PARALLEL_REPS else "")
          + (f"；失败重试 {getattr(a, 'retry', 0)} 次" if getattr(a, "retry", 0) else ""))
    gpu_sem = threading.Semaphore(gpu_jobs)
    t0 = time.time()
    results = []
    pending = list(todo)
    retries = int(getattr(a, "retry", 0) or 0)
    retry_wait = float(getattr(a, "retry_wait", 120) or 0)
    attempt = 0
    while True:
        attempt += 1
        if attempt > 1:
            print(f"\n第 {attempt} 轮（重试）：{len(pending)} 个作业", flush=True)
        wall0, mono0 = time.time(), time.monotonic()
        got = []
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = [pool.submit(worker, cfg, n, gpu_sem, None) for n in pending]
            for f in as_completed(futs):
                got.append(f.result())
        results.extend(got)
        # 机器睡过去时 time.time() 走得比 monotonic() 快；GROMACS 的 CUDA 调用
        # 常在恢复后报 "launch timed out"，所以把这层诊断打出来，便于判断重试
        # 是"偶发"还是"环境问题"。
        slept = (time.time() - wall0) - (time.monotonic() - mono0)
        if slept > 60:
            print(f"本轮实际耗时比单调时钟多 {slept / 60:.0f} 分钟：机器中途休眠/挂起过；"
                  f"GROMACS 可能在恢复时退出，重试会从 md.cpt 续跑", flush=True)
        bad_now = [r for r in got if not r[1]]
        if not bad_now or STOP.is_set() or attempt > retries:
            break
        pending = [r[0] for r in bad_now]
        wait = retry_wait * attempt
        if wait > 0:
            print(f"{wait:.0f} 秒后重试 {len(pending)} 个失败作业"
                  f"（第 {attempt}/{retries} 次重试）…", flush=True)
            end = time.monotonic() + wait
            while time.monotonic() < end and not STOP.is_set():
                time.sleep(min(5.0, max(0.0, end - time.monotonic())))
        if STOP.is_set():
            break
    print("队列执行结束")

    # 同一作业可能重试过多次，按最后一次的结果统计
    final = {r[0]: r for r in results}
    ok = [r for r in final.values() if r[1]]
    bad = [r for r in final.values() if not r[1]]
    hours = (time.time() - t0) / 3600.0
    summary = f"{len(ok)} 个完成" + (f"，{len(bad)} 个失败" if bad else "") \
              + f"；用时 {hours:.1f} 小时"
    if cfg.get("notify", True):
        if bad:
            names = "、".join(short_label(r[0], 18) for r in bad[:3])
            notify("MD 队列结束（有失败）", f"{summary}。失败：{names}", "fail")
        else:
            notify("MD 队列全部完成", summary, "info")


def cmd_status(a):
    q = queue()
    rows = []
    for n in q["jobs"]:
        meta = job_meta(n)
        if not meta:
            continue
        steps = [s for s in PREP_STEPS if done(n, s)]
        chk = meta.get("structure_check") or {}
        tag = chk.get("verdict", "-")
        if chk.get("repaired"):
            tag += "·已补"
        sc = job_verdict(job_dir(n)) or "-"
        rows.append((n, meta.get("status", "?"), len(steps), meta.get("prod_ns"),
                     meta.get("metrics", {}).get("em_max_force"), tag, sc))
    if not rows:
        print("队列为空")
        return
    print(f"{'作业':<24}{'状态':<10}{'准备':<6}{'生产(ns)':<10}{'体检':<8}{'自检':<6}{'EM最大力'}")
    for r in rows:
        print(f"{r[0]:<24}{str(r[1]):<10}{str(r[2]) + '/8':<6}"
              f"{str(r[3] or '-'):<10}{str(r[5]):<8}{r[6]:<6}"
              f"{r[4] if r[4] is not None else '-'}")


def cmd_selfcheck(a):
    q = queue()
    names = a.job or q["jobs"]
    bad = 0
    for n in names:
        meta = job_meta(n)
        if not meta:
            continue
        reps = int(meta.get("replicates", 1))
        for rep in range(1, reps + 1):
            sfx = f"_r{rep}" if reps > 1 else ""
            ana = job_dir(n) / f"07_analysis{sfx}"
            if not ana.exists():
                continue
            sc = check_replicate(ana, meta.get("prod_ns"))
            (ana / "selfcheck.json").write_text(
                json.dumps(sc, ensure_ascii=False, indent=2), encoding="utf-8")
            mark = "" if sc["verdict"] == "ok" else "  ←"
            print(f"[{n}] 重复{rep}: {sc['verdict']}{mark}")
            for c in sc["checks"]:
                if c["verdict"] != "ok":
                    print(f"    [{c['verdict']}] {c['check']}：{c['detail']}")
            bad += sc["verdict"] == "fail"
    print(f"自检完成{'，' + str(bad) + ' 个失败' if bad else ''}")


def cmd_report(a):
    import csv
    q = queue()
    out = PIPE / "summary.csv"
    keys = ["name", "status", "em_max_force", "nvt_T", "npt_density",
            "ns_per_day", "prod_ns", "started", "finished", "error"]
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for n in q["jobs"]:
            meta = job_meta(n)
            if not meta:
                continue
            row = {k: meta.get(k, "") for k in keys if k != "name"}
            row["name"] = n
            for k, v in meta.get("metrics", {}).items():
                row[k] = v
            w.writerow({k: row.get(k, "") for k in keys})
    print(f"已写出 {out}")


def cmd_remove(a):
    name = a.name
    q = queue()
    if name in q["jobs"]:
        q["jobs"].remove(name)
        save_queue(q)
    d = job_dir(name)
    if d.exists() and not a.keep:
        if input(f"确认删除 {d}? [y/N] ").lower() == "y":
            shutil.rmtree(d)
    print("已移除" if name not in q["jobs"] else "未找到")


# ---------------------------------------------------------------- 归档
# 每个作业里真正不可重建的只有 06_prod*/md.xtc 和拓扑/参数文件；下面这些
# 都是「占地方但能再算出来」的副本。默认只列出不删，加 --delete 才动手。
ARCHIVE_KINDS = [
    ("分析轨迹", ["07_analysis*/md_center.xtc", "07_analysis*/md_prot.xtc"],
     "mdpipe.py analysis --job {name}"),
    ("平衡轨迹", ["04_nvt*/nvt.xtc", "05_npt*/npt.xtc"],
     "平衡阶段已经结束，日常分析只用生产轨迹"),
]


def _sha256(path: Path, chunk: int = 4 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _size_h(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.2f} {unit}"
        n /= 1024.0
    return f"{n:.2f} TB"


def _archive_blockers(d: Path) -> dict[Path, str]:
    """返回 {文件: 不能删的原因}。只有当上游文件确实存在时才允许删副本。"""
    out: dict[Path, str] = {}
    for ana in sorted(d.glob("07_analysis*")):
        sfx = ana.name[len("07_analysis"):]          # "_r1" 或 ""
        prod = d / f"06_prod{sfx}"
        for fn in ("md_center.xtc", "md_prot.xtc"):
            p = ana / fn
            if not p.exists():
                continue
            if not (prod / "md.xtc").exists():
                out[p] = f"缺 {prod.name}/md.xtc，这份是唯一轨迹"
            elif not (prod / "md.tpr").exists():
                out[p] = f"缺 {prod.name}/md.tpr，重算分析需要它"
    for rep_dir in sorted(list(d.glob("04_nvt*")) + list(d.glob("05_npt*"))):
        xtc = rep_dir / ("nvt.xtc" if rep_dir.name.startswith("04_nvt") else "npt.xtc")
        if not xtc.exists():
            continue
        if rep_dir.name.startswith("04_nvt"):
            sfx = rep_dir.name[len("04_nvt"):]        # "_r1" 或 ""
            ok = (d / f"05_npt{sfx}" / "npt.gro").exists()
            why = f"缺 05_npt{sfx}/npt.gro，平衡还没做完"
        else:
            sfx = rep_dir.name[len("05_npt"):]
            ok = (d / f"06_prod{sfx}" / "md.xtc").exists()
            why = f"缺 06_prod{sfx}/md.xtc，生产还没开始"
        if not ok:
            out[xtc] = why
    return out


def _archive_plan(names, only_kinds=None, base: Path | None = None):
    """收集可删文件；顺带做完整性检查。"""
    plan = []
    for n in names:
        d = job_dir(n, base)
        if not d.is_dir():
            continue
        blocked = _archive_blockers(d)
        reps = len(list(d.glob("06_prod*"))) or 1
        for kind, patterns, note in ARCHIVE_KINDS:
            if only_kinds and kind not in only_kinds:
                continue
            for pat in patterns:
                for p in sorted(d.glob(pat)):
                    if p in blocked:
                        plan.append({"job": n, "kind": kind, "path": p,
                                     "size": p.stat().st_size,
                                     "blocked": blocked[p]})
                    else:
                        plan.append({"job": n, "kind": kind, "path": p,
                                     "size": p.stat().st_size, "blocked": None,
                                     "how": note.format(name=n, reps=reps)})
    return plan


def cmd_archive(a):
    base = Path(a.root).expanduser().resolve() if a.root else None
    if base:
        if not base.is_dir():
            sys.exit(f"作业根目录不存在：{base}")
        names = a.job or sorted(p.name for p in base.iterdir() if p.is_dir())
    else:
        q = queue()
        names = a.job or (q["jobs"] if a.all else [])
    if not names:
        sys.exit("需要 --job 作业名 或 --all")
    kinds = ["分析轨迹"] if a.analysis_only else (["平衡轨迹"] if a.eq_only else None)
    plan = _archive_plan(names, kinds, base)
    todo = [e for e in plan if not e["blocked"]]
    skip = [e for e in plan if e["blocked"]]
    if not plan:
        print("没有可归档的文件（分析轨迹和平衡轨迹都不在）")
        return

    # 目标副本校验：逐个比大小，--verify-hash 再比 sha256
    dest = Path(a.copy_to).expanduser().resolve() if a.copy_to else None
    bad_copy = []
    if dest:
        print(f"核对副本：{dest}")
        for e in todo:
            rel = e["path"].relative_to(job_dir(e["job"], base))
            tgt = dest / e["job"] / rel
            if not tgt.exists():
                bad_copy.append((e, f"目标缺 {tgt}"))
            elif tgt.stat().st_size != e["size"]:
                bad_copy.append((e, f"大小不一致 {tgt.stat().st_size} ≠ {e['size']}"))
            elif a.verify_hash:
                if sys.stdout.isatty():
                    print(f"  校验 {rel} …", end="\r", flush=True)
                if _sha256(tgt) != _sha256(e["path"]):
                    bad_copy.append((e, f"sha256 不一致 {tgt}"))
        if sys.stdout.isatty():
            print(" " * 60, end="\r")
        if bad_copy:
            print("副本校验未通过，以下文件不删：")
            for e, why in bad_copy:
                print(f"  {e['job']}/{e['path'].name}：{why}")
            ok = {id(e) for e, _ in bad_copy}
            skip += [{"job": e["job"], "path": e["path"], "kind": e["kind"],
                      "size": e["size"], "blocked": why} for e, why in bad_copy]
            todo = [e for e in todo if id(e) not in ok]
        else:
            print(f"  全部 {len(todo)} 个文件在目标位置核对通过")

    print(f"\n{'作业':<40}{'类别':<10}{'文件数':<8}{'体积'}")
    agg = {}
    for e in todo:
        k = (e["job"], e["kind"])
        c, s = agg.get(k, (0, 0))
        agg[k] = (c + 1, s + e["size"])
    for (n, kind), (c, s) in sorted(agg.items()):
        print(f"{n[:38]:<40}{kind:<10}{c:<8}{_size_h(s)}")
    total = sum(e["size"] for e in todo)
    print(f"{'合计':<40}{'':<10}{len(todo):<8}{_size_h(total)}")
    if skip:
        print(f"\n保留 {len(skip)} 个（有风险）：")
        for e in skip:
            print(f"  {e['job']}/{e['path'].name}：{e['blocked']}")

    if not a.delete:
        print("\n这是预演。确认无误后加 --delete 真正删除。")
        if dest:
            print("(--copy-to 只校验不复制；要复制先用系统命令把作业目录拷过去)")
        return
    if not todo:
        return

    deleted = 0
    for e in todo:
        try:
            e["path"].unlink()
            deleted += 1
        except OSError as exc:
            print(f"  删除失败 {e['path']}：{exc}")
    for d in {e["path"].parent for e in todo}:
        if d.name.startswith("07_analysis") and not any(d.iterdir()):
            d.rmdir()

    manifest = {
        "time": time.strftime("%Y-%m-%d %H:%M:%S"),
        "jobs": sorted({e["job"] for e in todo}),
        "copy_to": str(dest) if dest else None,
        "verified": bool(dest) and not bad_copy,
        "freed_bytes": sum(e["size"] for e in todo),
        "files": [{"job": e["job"], "kind": e["kind"],
                   "path": str(e["path"].relative_to(job_dir(e["job"], base))),
                   "bytes": e["size"]} for e in todo],
        "regenerate": {k: note for k, _, note in ARCHIVE_KINDS},
    }
    out = Path(a.manifest) if a.manifest else PIPE / (
        "archive_" + time.strftime("%Y%m%d_%H%M%S") + ".json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n已删除 {deleted} 个文件，释放 {_size_h(manifest['freed_bytes'])}")
    print(f"清单：{out}")


def short_label(name: str, limit: int = 26) -> str:
    """图例用的短名字。流水线作业名常是「课题前缀__配体_蛋白」，取后半段；
    中文长名会把图例撑到面板上，所以再按长度截断。"""
    lab = name.split("__")[-1].replace("_", " ")
    return lab if len(lab) <= limit else lab[:limit - 1] + "…"


def cmd_figures(a):
    """出图：单个作业出六联图，多个作业出对比图。"""
    sys.path.insert(0, str(ROOT / "scripts"))
    import make_figures

    names = a.job or (queue()["jobs"] if a.all else [])
    if not names:
        sys.exit("需要 --job 作业名 或 --all")
    names = [n for n in names if job_dir(n).exists()]
    if not names:
        sys.exit("作业目录不存在")
    if a.all:
        for n in names:
            make_figures.plot_job(job_dir(n), n, not a.no_pdf)
        return
    if len(names) == 1 and not a.compare:
        make_figures.plot_job(job_dir(names[0]), a.label or names[0], not a.no_pdf)
        return
    labels = a.labels or [short_label(n) for n in names]
    if len(labels) != len(names):
        sys.exit(f"--labels 给了 {len(labels)} 个，作业有 {len(names)} 个")
    outdir = Path(a.outdir) if a.outdir else PIPE / "figures"
    make_figures.plot_compare([job_dir(n) for n in names], labels, outdir, not a.no_pdf)


def cmd_viz(a):
    """复合物可视化：复合物 RMSD / 埋藏 SASA / 质心距离 + PyMOL 结构图。"""
    names = a.job or (queue()["jobs"] if a.all else [])
    if not names:
        sys.exit("需要 --job 作业名 或 --all")
    import viz
    cfg = config()
    for n in names:
        meta = job_meta(n)
        if not job_dir(n).exists() or not meta:
            print(f"跳过 {n}：没有作业目录或 job.json")
            continue
        viz.run_viz(cfg, job_dir(n), n, meta,
                    reps=[int(r) for r in a.replicates] if a.replicates else None,
                    frames=a.frames, want_pdf=not a.no_pdf,
                    do_pymol=not a.no_pymol)


def cmd_mmpbsa(a):
    """MM-PBSA 批量计算与对比图（转调 mmpbsa.py）。"""
    names = a.job or (queue()["jobs"] if a.all else [])
    if not names:
        sys.exit("需要 --job 作业名 或 --all")
    import mmpbsa
    argv = ["--job", *names, "--frames", str(a.frames),
            "--last-ns", str(a.last_ns)]
    if a.replicates:
        argv += ["--replicates", *map(str, a.replicates)]
    if a.outdir:
        argv += ["--outdir", a.outdir]
    if a.no_pdf:
        argv.append("--no-pdf")
    if a.keep_center:
        argv.append("--keep-center")
    mmpbsa.main(argv)


def cmd_interactions(a):
    """接触频率、氢键配对、配体聚类、2D 示意图（LigPlot+ 可选）。"""
    names = a.job or (queue()["jobs"] if a.all else [])
    if not names:
        sys.exit("需要 --job 作业名 或 --all")
    import interactions
    cfg = config()
    for n in names:
        meta = job_meta(n)
        if not job_dir(n).exists() or not meta:
            print(f"跳过 {n}：没有作业目录或 job.json")
            continue
        nrep = int(meta.get("replicates", 1))
        reps = [int(r) for r in a.replicates] if a.replicates else list(range(1, nrep + 1))
        for rep in reps:
            print(f"[{n}] 相互作用 重复 {rep} …", flush=True)
            interactions.run_job(
                cfg, n, rep, meta, cutoff=a.cutoff, dt_ns=a.dt,
                want_pdf=not a.no_pdf, do_map=not a.no_map,
                do_ligplot=not a.no_ligplot)


def cmd_analysis(a):
    """按当前代码重算分析：轨迹修正 + 六项指标 + PCA + 出图，不重跑 MD。

    分析代码改过之后，早先跑完的作业还留着旧口径的 xvg/图片。这条命令用磁盘上的
    新代码把这些作业重新处理一遍——轨迹（06_prod*/md.xtc）不动，只重算分析。
    """
    names = a.job or (queue()["jobs"] if a.all else [])
    if not names:
        sys.exit("需要 --job 作业名 或 --all")
    names = [n for n in names if job_dir(n).exists()]
    if not names:
        sys.exit("作业目录不存在")

    cfg = config()
    for n in names:
        meta = job_meta(n)
        if not meta:
            print(f"[{n}] 跳过：没有 job.json", flush=True)
            continue
        nrep = meta.get("replicates", 1)
        reps = [int(r) for r in a.replicates] if a.replicates else list(range(1, nrep + 1))
        logs = job_dir(n) / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        any_done = False
        for rep in reps:
            sfx = f"_r{rep}" if nrep > 1 else ""
            if not (job_dir(n) / f"06_prod{sfx}" / "md.xtc").exists():
                print(f"[{n}] 跳过 r{rep}：没有 06_prod{sfx}/md.xtc", flush=True)
                continue
            analysis_step(cfg, n, meta, logs / f"analysis_r{rep}.log", rep)
            pca_step(cfg, n, meta, logs / f"pca_r{rep}.log", rep, force=a.force)
            any_done = True
        if any_done:
            figures_step(cfg, n, meta, logs / "figures.log")
            # 顺手清掉旧口径留下的轨迹，免得以后又拿错文件
            for stale in job_dir(n).glob("07_analysis*/md_nojump.xtc"):
                stale.unlink(missing_ok=True)
                print(f"[{n}] 已删除旧轨迹 {stale.relative_to(job_dir(n))}", flush=True)


def cmd_fel(a):
    """PCA + 自由能景观：先算（必要时强制），再出图。"""
    names = a.job or (queue()["jobs"] if a.all else [])
    if not names:
        sys.exit("需要 --job 作业名 或 --all")
    names = [n for n in names if job_dir(n).exists()]
    if not names:
        sys.exit("作业目录不存在")

    cfg = config()
    logs = job_dir(names[0]) / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    for n in names:
        meta = job_meta(n)
        reps = list(range(1, meta.get("replicates", 1) + 1))
        for rep in reps:
            sfx = f"_r{rep}" if len(reps) > 1 else ""
            if not a.force and (job_dir(n) / f"07_analysis{sfx}" / "pca" / "2dproj.xvg").exists():
                continue
            pca_step(cfg, n, meta, job_dir(n) / "logs" / f"pca_r{rep}.log", rep,
                     force=a.force)

    sys.path.insert(0, str(ROOT / "scripts"))
    import make_fel
    if a.compare and len(names) > 1:
        outdir = Path(a.outdir) if a.outdir else PIPE / "figures"
        labels = a.labels or [short_label(n) for n in names]
        make_fel.plot_compare([job_dir(n) for n in names], labels, outdir, not a.no_pdf)
        return
    for n in names:
        make_fel.plot_job(job_dir(n), n, not a.no_pdf)


def cmd_gui(a):
    """从任意目录启动图形界面：python <本文件全路径> gui"""
    sys.path.insert(0, str(ROOT / "scripts"))
    import mdpipe_gui
    mdpipe_gui.main()


def cmd_doctor(a):
    """环境自检；返回码非 0 表示有必需项缺失。"""
    sys.path.insert(0, str(ROOT / "scripts"))
    import doctor
    argv = []
    for flag in ("plan", "json", "install_pip", "with_optional", "yes"):
        if getattr(a, flag, False):
            argv.append("--" + flag.replace("_", "-"))
    sys.exit(doctor.main(argv))


CONFIG_KEYS = ("gmx", "wsl_distro", "wsl_home", "amber_home", "mmpbsa_env",
               "gmx_share", "ligplot_home")


def cmd_config(a):
    """查看或探测组件路径；安装程序用 --detect --write 写首次配置。"""
    import autoconfig
    cfg = config()
    if not a.detect:
        for k in CONFIG_KEYS:
            print(f"{k:<14}{cfg.get(k, '') or '（未设置）'}")
        print(f"\n配置文件：{CONFIG_PATH}")
        print("探测本机路径：mdpipe config --detect [--write]")
        return

    det = autoconfig.detect_all(cfg)
    print(f"{'键':<14}{'当前值':<42}{'探测结果'}")
    print("-" * 92)
    for k in CONFIG_KEYS:
        cur = str(cfg.get(k, "")) or "（未设置）"
        new = det.get(k, "")
        print(f"{k:<14}{cur[:40]:<42}{new or '（未探测到）'}")

    if not a.write:
        print("\n只显示探测结果。加 --write 写入 config.json。")
        return
    merged = autoconfig.merge_detected(cfg, det, overwrite=a.overwrite)
    save_json(CONFIG_PATH, merged)
    changed = [k for k in det if merged.get(k) != cfg.get(k)]
    print(f"\n已写入 {CONFIG_PATH}")
    print("更新：" + ("、".join(changed) if changed else "无变化"))


DOC_FILES = {
    "manual": "分子动力学操作手册.md",
    "checklist": "安装与环境清单.md",
    "install": "安装说明.md",
    "readme": "README.md",
}


def doc_path(key: str) -> Path:
    """文档在安装版是 {app}\\docs，源码版散在根目录 / installer / pipeline。"""
    name = DOC_FILES[key]
    if key == "install":
        cands = [ROOT / "docs" / name, ROOT / "installer" / name]
    elif key == "readme":
        cands = [ROOT / "pipeline" / name, ROOT / "docs" / name]
    else:
        cands = [ROOT / "docs" / name, ROOT / name]
    for p in cands:
        if p.is_file():
            return p
    return cands[0]


def _doc_section(lines: list[str], sec: str) -> list[str]:
    """取某一节：支持 16 或 20.10 这种编号。"""
    pat = re.compile(r"^(#{2,4})\s+" + re.escape(sec) + r"[.．\s]")
    start = level = None
    for i, line in enumerate(lines):
        m = pat.match(line)
        if m:
            start, level = i, len(m.group(1))
            break
    if start is None:
        return []
    out = [lines[start]]
    for line in lines[start + 1:]:
        m = re.match(r"^(#{2,4})\s", line)
        if m and len(m.group(1)) <= level:
            break
        out.append(line)
    return out


def cmd_doc(a):
    """列出/检索手册与清单，供人和 AI 查完整说明。"""
    if a.list or not a.name:
        print(f"{'名字':<10}{'状态':<6}路径")
        for k in DOC_FILES:
            p = doc_path(k)
            print(f"{k:<10}{'有' if p.is_file() else '缺':<6}{p}")
        print("\n用法：mdpipe doc manual --section 16 | --grep GROMOS")
        return

    p = doc_path(a.name)
    if not p.is_file():
        sys.exit(f"找不到文档：{p}")
    text = p.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    if a.grep:
        pat = re.compile(a.grep, re.I)
        hits = [i for i, line in enumerate(lines) if pat.search(line)]
        if not hits:
            print(f"{p} 里没有匹配 {a.grep!r} 的行")
            return
        for i in hits:
            lo, hi = max(0, i - a.context), min(len(lines), i + a.context + 1)
            for j in range(lo, hi):
                mark = ">" if j == i else " "
                print(f"{mark}{j + 1:>6}: {lines[j]}")
            print()
        return

    if a.section:
        chunk = _doc_section(lines, a.section)
        if not chunk:
            sys.exit(f"{p.name} 里没有第 {a.section} 节")
        print("\n".join(chunk))
        return

    print(f"{p}\n（{len(lines)} 行，{len(text)} 字符）")
    print("用 --grep 关键词 或 --section 编号 取具体内容，别整份读进上下文。")


def cmd_skill(a):
    """把 md-pipeline 技能装到指定 AI Harness 的目录（转调 install_skill.py）。"""
    sys.path.insert(0, str(ROOT / "scripts"))
    import install_skill
    argv = []
    if a.target:
        argv += ["--target", *a.target]
    if a.list:
        argv.append("--list")
    if a.print:
        argv.append("--print")
    if a.force:
        argv.append("--force")
    sys.exit(install_skill.main(argv))


def main():
    ap = argparse.ArgumentParser(description="蛋白–小分子 MD 批处理流水线")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("add", help="登记作业（可一对多、多对多）")
    p.add_argument("--protein", help="蛋白 PDB")
    p.add_argument("--protein-dir", dest="protein_dir", help="蛋白目录：目录下每个 PDB 一个蛋白")
    p.add_argument("--ligand", nargs="+", help="配体位姿，可以给多个文件")
    p.add_argument("--ligand-dir", dest="ligand_dir", help="配体目录：目录下每个分子一个作业")
    p.add_argument("--complex", nargs="+", metavar="复合物.pdb",
                   help="直接投「蛋白+配体」复合物 PDB，自动拆成两份")
    p.add_argument("--lig-resname", dest="lig_resname",
                   help="复合物里配体的残基名；不给则自动判断")
    p.add_argument("--name", help="作业名；一对多时作为前缀，默认取文件名")
    p.add_argument("--charge", type=int, help="净电荷；不给则自动判定")
    p.add_argument("--resname", default="LIG")
    p.add_argument("--model", type=int, default=1, help="取第几个对接构象")
    p.add_argument("--replicates", type=int, default=1)
    p.add_argument("--ns", type=float, default=200, help="生产时长 ns")
    p.add_argument("--chain", help="只保留该链")
    p.add_argument("--force-altloc", dest="force_altloc", help="如 120:B,127:B")
    p.add_argument("--keep-het", dest="keep_het", help="保留的 HETATM 残基名，逗号分隔，如 ZN,MG")
    p.add_argument("--prepared", action="store_true", help="蛋白 PDB 已清理过，跳过预处理")
    p.add_argument("--seq-dir", dest="seq_dir",
                   help="FASTA 序列目录（{蛋白名}.fasta）。受体文件没有 SEQRES 时，"
                        "自动补链需要它；序列从 UniProt 下载")
    p.add_argument("--seqres-dir", dest="seqres_dir",
                   help="原版 PDB 目录，给 *_clean.pdb 借 SEQRES")
    p.add_argument("--forcefield", help="覆盖 config 的力场，如 oplsaa")
    p.add_argument("--water", help="覆盖 config 的水模型，如 spc")
    p.add_argument("--ligand-itp", dest="ligand_itp",
                   help="外部配体拓扑（LigParGen / ATB / CGenFF），给了就不再跑 antechamber")
    p.add_argument("--ligand-gro", dest="ligand_gro",
                   help="配体坐标 gro，配合 --ligand-itp 使用")
    p.add_argument("--ligand-posre", dest="ligand_posre",
                   help="配体的位置约束 itp；不给则按 gro 自动生成")
    p.add_argument("--ligand-ff", dest="ligand_ff",
                   help="配体参数的来源（gaff2/opls/cgenff），用于力场混用检查")
    p.add_argument("--ligand-atb", dest="ligand_atb",
                   help="ATB 导出的配体参数包（目录或 zip/tgz），自动整理后替代 antechamber")
    p.add_argument("--pair", action="append", metavar="蛋白::配体1,配体2",
                   help='指定配对，可写多次；例 --pair "FABP5::Glc,Gal" --pair "FABP6::Man"')
    p.add_argument("--csv", help="配对表 CSV：列 protein,ligand[,name,charge,replicates,ns]")
    p.add_argument("--yes", action="store_true", help="超过 50 个作业时确认")
    p.set_defaults(func=cmd_add)

    p = sub.add_parser("complex", help="复合物 PDB 专用入口（自动拆成蛋白 + 配体）")
    p.add_argument("files", nargs="+", metavar="复合物.pdb",
                   help="一个或多个「蛋白 + 配体」复合物 PDB")
    p.add_argument("--name", help="作业名前缀；多个复合物时自动拼上各自文件名")
    p.add_argument("--lig-resname", dest="lig_resname",
                   help="配体的残基名；不给则自动判断")
    p.add_argument("--charge", type=int, help="净电荷；不给则自动判定")
    p.add_argument("--resname", default="LIG")
    p.add_argument("--model", type=int, default=1, help="取第几个对接构象")
    p.add_argument("--replicates", type=int, default=1)
    p.add_argument("--ns", type=float, default=None, help="生产时长 ns（默认取 config）")
    p.add_argument("--chain", help="只保留该链")
    p.add_argument("--force-altloc", dest="force_altloc", help="如 120:B,127:B")
    p.add_argument("--keep-het", dest="keep_het", help="保留的 HETATM 残基名，逗号分隔")
    p.add_argument("--prepared", action="store_true", help="蛋白部分已清理过，跳过预处理")
    p.add_argument("--seq-dir", dest="seq_dir",
                   help="FASTA 序列目录（{蛋白名}.fasta）。受体文件没有 SEQRES 时，"
                        "自动补链需要它；序列从 UniProt 下载")
    p.add_argument("--seqres-dir", dest="seqres_dir",
                   help="原版 PDB 目录，给 *_clean.pdb 借 SEQRES")
    p.add_argument("--forcefield", help="覆盖 config 的力场，如 oplsaa")
    p.add_argument("--water", help="覆盖 config 的水模型，如 spc")
    p.add_argument("--ligand-itp", dest="ligand_itp",
                   help="外部配体拓扑（LigParGen / ATB / CGenFF），给了就不再跑 antechamber")
    p.add_argument("--ligand-gro", dest="ligand_gro",
                   help="配体坐标 gro，配合 --ligand-itp 使用")
    p.add_argument("--ligand-posre", dest="ligand_posre",
                   help="配体的位置约束 itp；不给则按 gro 自动生成")
    p.add_argument("--ligand-ff", dest="ligand_ff",
                   help="配体参数的来源（gaff2/opls/cgenff），用于力场混用检查")
    p.add_argument("--ligand-atb", dest="ligand_atb",
                   help="ATB 导出的配体参数包（目录或 zip/tgz），自动整理后替代 antechamber")
    p.add_argument("--yes", action="store_true", help="超过 50 个作业时确认")
    p.set_defaults(func=cmd_complex)

    p = sub.add_parser("scan", help="扫描 pipeline/inbox/ 自动登记")
    p.set_defaults(func=cmd_scan)

    p = sub.add_parser("run", help="执行队列")
    p.add_argument("--workers", type=int, default=None, help="同时处理的作业数（默认 2）")
    p.add_argument("--gpu-jobs", dest="gpu_jobs", type=int, default=None,
                   help="同时占用 GPU 的 mdrun 数（默认 2，本机实测与单个跑总吞吐持平）")
    p.add_argument("--job", nargs="*", help="只跑指定作业")
    p.add_argument("--ns", type=float, help="覆盖生产时长")
    p.add_argument("--nvt-ps", dest="nvt_ps", type=float, help="覆盖 NVT 时长")
    p.add_argument("--npt-ps", dest="npt_ps", type=float, help="覆盖 NPT 时长")
    p.add_argument("--force", action="store_true", help="已完成的作业也重跑")
    p.add_argument("--prep-only", dest="prep_only", action="store_true",
                   help="只做结构准备（拓扑/溶剂/离子/索引），不占 GPU")
    p.add_argument("--only-failed", dest="only_failed", action="store_true")
    p.add_argument("--retry", type=int, default=0,
                   help="失败作业自动重试的次数（默认 0）；重试会从 md.cpt 续跑，"
                        "适合笔记本中途休眠后 CUDA 调用超时的情形")
    p.add_argument("--retry-wait", dest="retry_wait", type=float, default=120,
                   help="首次重试前等待的秒数，之后逐轮翻倍（默认 120）")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--parallel-replicates", dest="parallel_replicates", action="store_true",
                   default=None, help="同一作业的多个重复并行跑（默认在作业数少于并行槽位时自动启用）")
    g.add_argument("--serial-replicates", dest="parallel_replicates", action="store_false",
                   help="强制重复顺序跑")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("status", help="查看队列状态")
    p.set_defaults(func=cmd_status)

    p = sub.add_parser("selfcheck", help="分析结果自检（只读 xvg，不重算）")
    p.add_argument("--job", nargs="+", help="作业名，可给多个")
    p.add_argument("--all", action="store_true", help="队列里所有作业")
    p.set_defaults(func=cmd_selfcheck)

    p = sub.add_parser("archive", help="清理可重算轨迹（默认只列出，--delete 才删）")
    p.add_argument("--job", nargs="+", help="作业名，可给多个")
    p.add_argument("--root", help="作业根目录，默认 pipeline/jobs；指向别的结果包时用")
    p.add_argument("--all", action="store_true", help="队列里所有作业")
    p.add_argument("--delete", action="store_true", help="真正删除；不加只做预演")
    p.add_argument("--copy-to", dest="copy_to",
                   help="完整结果包的目录（如 E:\\分子动力学模拟）；只校验副本，不复制")
    p.add_argument("--verify-hash", dest="verify_hash", action="store_true",
                   help="校验 sha256 而不只是文件大小（大文件会慢）")
    p.add_argument("--analysis-only", dest="analysis_only", action="store_true",
                   help="只归档分析轨迹（md_center/md_prot.xtc）")
    p.add_argument("--eq-only", dest="eq_only", action="store_true",
                   help="只归档平衡轨迹（nvt/npt.xtc）")
    p.add_argument("--manifest", help="清单输出路径，默认 pipeline/archive_<时间>.json")
    p.set_defaults(func=cmd_archive)

    p = sub.add_parser("provenance", help="看作业的参数溯源（gmx/AmberTools/mdp 哈希）")
    p.add_argument("--job", nargs="+", help="作业名，可给多个")
    p.add_argument("--all", action="store_true", help="队列里所有作业")
    p.add_argument("--backfill", action="store_true",
                   help="给没有溯源信息的老作业补录（标记 backfilled）")
    p.set_defaults(func=cmd_provenance)

    p = sub.add_parser("report", help="汇总为 pipeline/summary.csv")
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("remove", help="移除作业")
    p.add_argument("name")
    p.add_argument("--keep", action="store_true", help="只移出队列，保留目录")
    p.set_defaults(func=cmd_remove)

    p = sub.add_parser("fel", help="PCA 自由能景观（计算 + 出图）")
    p.add_argument("--job", nargs="+", help="作业名，可给多个")
    p.add_argument("--all", action="store_true", help="队列里所有作业")
    p.add_argument("--compare", action="store_true", help="多个作业画进同一张 3D 对比图")
    p.add_argument("--force", action="store_true", help="轨迹短于 fel_min_ns 也照算")
    p.add_argument("--labels", nargs="+", help="对比图图例，按作业顺序给")
    p.add_argument("--outdir", help="对比图输出目录，默认 pipeline/figures")
    p.add_argument("--no-pdf", action="store_true", help="只出 PNG")
    p.set_defaults(func=cmd_fel)

    p = sub.add_parser("figures", help="出图（单体系六联图 / 多体系对比图）")
    p.add_argument("--job", nargs="+", help="作业名，可给多个")
    p.add_argument("--all", action="store_true", help="队列里所有作业各出一张")
    p.add_argument("--compare", action="store_true", help="多个作业画进同一张对比图")
    p.add_argument("--label", help="单体系图里的显示名")
    p.add_argument("--labels", nargs="+", help="对比图图例，按作业顺序给")
    p.add_argument("--outdir", help="对比图输出目录，默认 pipeline/figures")
    p.add_argument("--no-pdf", action="store_true", help="只出 PNG")
    p.set_defaults(func=cmd_figures)

    p = sub.add_parser("viz", help="复合物可视化：补充指标 + PyMOL 快照/B-factor")
    p.add_argument("--job", nargs="+", help="作业名，可给多个")
    p.add_argument("--all", action="store_true", help="队列里所有作业")
    p.add_argument("--replicates", nargs="+", help="只做指定重复，如 1 2")
    p.add_argument("--frames", nargs="+", type=float,
                   help="快照时间 (ns)，默认 0 / 中 / 末")
    p.add_argument("--no-pymol", action="store_true", help="只算三条曲线，不渲染结构")
    p.add_argument("--no-pdf", action="store_true", help="只出 PNG")
    p.set_defaults(func=cmd_viz)

    p = sub.add_parser("interactions",
                       help="配体-蛋白相互作用：接触/氢键/聚类/2D 示意图")
    p.add_argument("--job", nargs="+", help="作业名，可给多个")
    p.add_argument("--all", action="store_true", help="队列里所有作业")
    p.add_argument("--replicates", nargs="+", help="只做指定重复，如 1 2")
    p.add_argument("--cutoff", type=float, default=0.45,
                   help="接触截断，单位 nm（默认 0.45）")
    p.add_argument("--dt", type=float, default=0.2, help="采样间隔 ns（默认 0.2）")
    p.add_argument("--no-map", action="store_true", help="不画 RDKit 2D 示意图")
    p.add_argument("--no-ligplot", action="store_true",
                   help="不调用 LigPlot+（装了就默认调用）")
    p.add_argument("--no-pdf", action="store_true", help="只出 PNG")
    p.set_defaults(func=cmd_interactions)

    p = sub.add_parser("mmpbsa", help="MM-PBSA 结合自由能（批量 + 对比图）")
    p.add_argument("--job", nargs="+", help="作业名，可给多个")
    p.add_argument("--all", action="store_true", help="队列里所有作业")
    p.add_argument("--replicates", nargs="+", help="只算指定重复，如 1 2")
    p.add_argument("--frames", type=int, default=10, help="抽多少帧（默认 10）")
    p.add_argument("--last-ns", type=float, default=20.0,
                   help="取样窗口：轨迹最后 N ns（默认 20，0=全程）")
    p.add_argument("--keep-center", action="store_true",
                   help="保留居中轨迹 center.xtc（默认算完删除）")
    p.add_argument("--outdir", help="对比图目录，默认 pipeline/figures")
    p.add_argument("--no-pdf", action="store_true", help="只出 PNG")
    p.set_defaults(func=cmd_mmpbsa)

    p = sub.add_parser("analysis", help="按当前代码重算分析与图（不重跑 MD）")
    p.add_argument("--job", nargs="+", help="作业名，可给多个")
    p.add_argument("--all", action="store_true", help="队列里所有作业")
    p.add_argument("--replicates", nargs="+", help="只重算指定重复，如 1 2")
    p.add_argument("--force", action="store_true", help="轨迹短于 fel_min_ns 也照算 PCA")
    p.set_defaults(func=cmd_analysis)

    p = sub.add_parser("config", help="查看/探测组件路径（安装时自动写入）")
    p.add_argument("--detect", action="store_true",
                   help="探测本机 GROMACS / WSL / AmberTools / MM-PBSA / LigPlot+")
    p.add_argument("--write", action="store_true", help="把探测结果写入 config.json")
    p.add_argument("--overwrite", action="store_true",
                   help="覆盖已有路径（默认只补空项）")
    p.set_defaults(func=cmd_config)

    p = sub.add_parser("doctor", help="环境自检：GROMACS / WSL / AmberTools / 出图工具")
    p.add_argument("--plan", action="store_true",
                   help="只打印安装清单（哪些能自动装、哪些必须手动），不检测")
    p.add_argument("--json", action="store_true", help="机器可读输出")
    p.add_argument("--install-pip", dest="install_pip", action="store_true",
                   help="自动安装缺失的必需 Python 包（numpy/scipy/matplotlib；需网络）")
    p.add_argument("--with-optional", dest="with_optional", action="store_true",
                   help="配合 --install-pip 一并装 pandas/rdkit/PyMOL")
    p.add_argument("--yes", action="store_true", help="安装前不再询问")
    p.set_defaults(func=cmd_doctor)

    p = sub.add_parser("doc", help="查看/检索操作手册、安装清单（给人和 AI 查细节）")
    p.add_argument("name", nargs="?", choices=list(DOC_FILES),
                   help="manual=操作手册，checklist=安装清单，install=安装说明，readme=流水线说明")
    p.add_argument("--list", action="store_true", help="列出四份文档的路径")
    p.add_argument("--grep", help="按关键词检索（忽略大小写），给出行号")
    p.add_argument("--section", help="打印某一节，如 16 或 20.10")
    p.add_argument("--context", type=int, default=2, help="grep 前后各显示几行，默认 2")
    p.set_defaults(func=cmd_doc)

    p = sub.add_parser("skill", help="把 AI 操作技能装到 ClawsGO / Claude / AGENTS.md")
    p.add_argument("--target", nargs="+",
                   choices=["clawsgo", "dsh", "workbuddy", "modao", "agentshome",
                            "codex", "claude", "project", "pclaude", "agents",
                            "trae", "qoder", "cursor", "copilot"],
                   help="安装目标，可给多个；detected=本机检测到的全部")
    p.add_argument("--list", action="store_true", help="列出各目标的状态")
    p.add_argument("--print", dest="print", action="store_true", help="打印技能内容")
    p.add_argument("--force", action="store_true", help="覆盖已存在的不同版本")
    p.set_defaults(func=cmd_skill)

    p = sub.add_parser("gui", help="打开图形界面")
    p.set_defaults(func=cmd_gui)

    a = ap.parse_args()
    if a.cmd == "add" and not (a.csv or a.pair or a.complex or a.protein or a.protein_dir):
        ap.error("需要 --complex（复合物 PDB）、--protein/--ligand，或用 --pair / --csv")
    a.func(a)


if __name__ == "__main__":
    main()
