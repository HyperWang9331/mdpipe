#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mdpipe_gui.py —— mdpipe 的图形界面（Tkinter，无需安装任何东西）

    python scripts/mdpipe_gui.py

左侧填文件、点"添加到队列"；右侧点"开始运行"。所有耗时操作都在后台线程，
界面不会卡；日志实时显示，"停止"会在当前这一小步结束后生效。
"""
from __future__ import annotations

import contextlib
import io
import os
import queue
import shutil
import subprocess
import sys
import threading
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from tkinter import (BooleanVar, DoubleVar, Listbox, Menu, StringVar, Text, Tk,
                     Toplevel, filedialog, messagebox)
from tkinter import ttk

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdpipe  # noqa: E402
import mdview  # noqa: E402

APP_VERSION = "1.7"


def _find_doc(name: str, root: Path | None = None) -> Path | None:
    """安装版文档在 {app}\\docs，源码版在仓库根目录，两处都要找。"""
    root = root or mdpipe.ROOT
    for base in (root / "docs", root, root / "pipeline", root / "installer",
                 root / "skills" / "md-pipeline"):
        p = base / name
        if p.is_file():
            return p
    return None


def _icon_path() -> Path | None:
    """找应用图标：安装版在 {app}\appicon.ico，源码版在 pipeline/figures/。"""
    here = Path(__file__).resolve().parent
    for p in (here.parent / "appicon.ico",
              here.parent / "pipeline" / "figures" / "appicon.ico",
              here / "appicon.ico"):
        if p.is_file():
            return p
    return None

# 可在「设置」里改的组件路径；顺序就是对话框里的顺序
PATH_FIELDS = [
    ("gmx", "GROMACS", "点「浏览…」选 gmx.exe，或点「自动探测」", True),
    ("wsl_distro", "WSL 发行版", "如 Ubuntu-24.04", False),
    ("amber_home", "AmberTools", "conda 环境路径，如 ~/miniconda3/envs/amber", False),
    ("mmpbsa_env", "gmx_MMPBSA", "conda 环境路径，如 ~/miniconda3/envs/mmpbsa", False),
    ("gmx_share", "GMXLIB", "留空则从 GROMACS 路径推算", False),
    ("wsl_home", "WSL 家目录", "如 /home/用户名", False),
    ("ligplot_home", "LigPlot+ 根目录", "留空则自动查找，或设环境变量 LIGPLUS_HOME", True),
]


ABOUT_TEXT = """MDpipe  v{version}
蛋白–配体分子动力学流水线（Windows + WSL）

把蛋白–小分子分子动力学模拟的结构准备、生产运行、分析和出图
串成一条可复现的流水线。

───────── 计算与能量 ─────────
GROMACS
    https://www.gromacs.org        LGPL-2.1
    Abraham et al., SoftwareX 1 (2015) 19–25

AmberTools / ParmEd / cpptraj
    https://ambermd.org

gmx_MMPBSA 1.7
    https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA
    Valdés-Tresanco et al., J. Chem. Theory Comput. 17 (2021) 6281–6291
    MM-PBSA 方法：Kollman et al., Acc. Chem. Res. 33 (2000) 889–897

OpenBabel
    https://openbabel.org           GPL-2.0

───────── 分析与绘图 ─────────
NumPy · SciPy · pandas            BSD-3
Matplotlib                        https://matplotlib.org
RDKit                             https://www.rdkit.org        BSD-3
PyMOL                             https://pymol.org
LigPlot+ 2.2
    https://www.ebi.ac.uk/thornton-srv/software/LigPlus/
    Laskowski & Swindells, J. Chem. Inf. Model. 51 (2011) 2778–2786
Ghostscript                       https://www.ghostscript.com   AGPL-3.0

───────── 力场与参数 ─────────
Amber14SB / TIP3P                 Case et al., Amber 2024
GROMOS 54A7 + ATB
    https://atb.uq.edu.au
    Malde et al., J. Comput. Chem. 32 (2011) 2497–2509

各组件按其自身许可证使用；本程序不重新分发这些组件的二进制。
安装后可用「设置 → 运行环境自检」确认每一项是否就位。

本工具负责把流程跑完并给出可核对的数字，结论仍需人来判断。
"""


class QueueWriter:
    """把 print 的输出送进界面日志。"""

    def __init__(self, q: queue.Queue):
        self.q = q

    def write(self, text):
        if text:
            self.q.put(text)

    def flush(self):
        pass


class PairDialog(Toplevel):
    """列出所有"蛋白 × 配体"组合，默认全选，点一行切换要不要跑。"""

    def __init__(self, master, pairs, title="选择要跑的配对"):
        super().__init__(master)
        self.title(title)
        self.geometry("640x480")
        self.pairs = list(pairs)
        self.chosen = set(range(len(self.pairs)))
        self.result: list | None = None

        ttk.Label(self, text="点一行切换勾选；默认全选。只跑需要的组合。",
                  padding=(10, 8)).pack(anchor="w")

        frm = ttk.Frame(self)
        frm.pack(fill="both", expand=True, padx=10)
        self.tree = ttk.Treeview(frm, columns=("sel", "prot", "lig"), show="headings")
        self.tree.heading("sel", text="跑")
        self.tree.heading("prot", text="蛋白")
        self.tree.heading("lig", text="配体")
        self.tree.column("sel", width=50, anchor="center")
        self.tree.column("prot", width=250)
        self.tree.column("lig", width=250)
        self.tree.pack(side="left", fill="both", expand=True)
        sb = ttk.Scrollbar(frm, orient="vertical", command=self.tree.yview)
        sb.pack(side="right", fill="y")
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.bind("<Button-1>", self._on_click)
        self._fill()

        bar = ttk.Frame(self)
        bar.pack(fill="x", padx=10, pady=10)
        ttk.Button(bar, text="全选", command=lambda: self._set_all(True)).pack(side="left")
        ttk.Button(bar, text="全不选", command=lambda: self._set_all(False)).pack(side="left", padx=4)
        ttk.Button(bar, text="反选", command=self._invert).pack(side="left")
        ttk.Label(bar, text="").pack(side="left", expand=True)
        ttk.Button(bar, text="取消", command=self.destroy).pack(side="right")
        ttk.Button(bar, text="确定", command=self._ok).pack(side="right", padx=4)

    def _fill(self):
        for i in self.tree.get_children():
            self.tree.delete(i)
        for i, (p, l) in enumerate(self.pairs):
            self.tree.insert("", "end", iid=str(i),
                             values=("☑" if i in self.chosen else "☐", p.name, l.name))

    def _refresh_row(self, i):
        sel, prot, lig = self.tree.item(str(i))["values"]
        self.tree.item(str(i), values=("☑" if i in self.chosen else "☐", prot, lig))

    def _on_click(self, event):
        row = self.tree.identify_row(event.y)
        if not row:
            return
        i = int(row)
        if i in self.chosen:
            self.chosen.discard(i)
        else:
            self.chosen.add(i)
        self._refresh_row(i)

    def _set_all(self, on: bool):
        self.chosen = set(range(len(self.pairs))) if on else set()
        self._fill()

    def _invert(self):
        self.chosen = set(range(len(self.pairs))) - self.chosen
        self._fill()

    def _ok(self):
        self.result = [self.pairs[i] for i in sorted(self.chosen)]
        self.destroy()

    def show(self):
        self.transient(self.master)
        self.grab_set()
        self.wait_window(self)
        return self.result


class App:
    def __init__(self, root: Tk):
        self.root = root
        self.cfg = mdpipe.config()
        self.logq: queue.Queue = queue.Queue()
        self.worker_thread: threading.Thread | None = None
        self.running = False
        self.ligand_list: list[str] = []

        root.title("分子动力学批处理 · MDpipe")
        root.geometry("980x680")
        root.minsize(860, 600)

        icon = _icon_path()
        if icon:
            try:
                root.iconbitmap(default=str(icon))
            except Exception:
                pass

        self._build_menu()
        self._build_form()
        self._build_table()
        self._build_log()

        self.root.after(200, self._drain_log)
        self.refresh()

    # ---------------- 界面搭建 ----------------
    def _build_form(self):
        frm = ttk.LabelFrame(self.root, text="新建作业")
        frm.pack(fill="x", padx=10, pady=(10, 6))

        self.var_protein = StringVar()
        self.var_complex = StringVar()
        self.var_ligand = StringVar()
        self.var_name = StringVar()
        self.var_charge = StringVar(value="0")
        self.var_ns = StringVar(value="200")
        self.var_reps = StringVar(value="1")
        self.var_workers = StringVar(value=str(self.cfg.get("workers", 2)))
        self.var_gpu = StringVar(value=str(self.cfg.get("gpu_jobs", 2)))
        self.var_auto_charge = BooleanVar(value=False)

        r1 = ttk.Frame(frm)
        r1.pack(fill="x", padx=8, pady=4)
        ttk.Label(r1, text="蛋白 PDB").pack(side="left")
        ttk.Entry(r1, textvariable=self.var_protein, width=52).pack(side="left", padx=6)
        ttk.Button(r1, text="浏览…", command=self._pick_protein).pack(side="left")
        ttk.Label(r1, text="（只放蛋白；蛋白+配体一起的走下面「复合物」）").pack(side="left", padx=6)

        rc = ttk.Frame(frm)
        rc.pack(fill="x", padx=8, pady=4)
        ttk.Label(rc, text="复合物 PDB").pack(side="left")
        ttk.Entry(rc, textvariable=self.var_complex, width=50).pack(side="left", padx=6)
        ttk.Button(rc, text="选文件…", command=self._pick_complex).pack(side="left")
        ttk.Button(rc, text="添加复合物", command=self.on_add_complex).pack(side="left", padx=(6, 0))
        ttk.Label(rc, text="（可多选，自动拆成蛋白+配体）").pack(side="left", padx=6)

        r2 = ttk.Frame(frm)
        r2.pack(fill="x", padx=8, pady=4)
        ttk.Label(r2, text="配体位姿").pack(side="left")
        ttk.Entry(r2, textvariable=self.var_ligand, width=44).pack(side="left", padx=6)
        ttk.Button(r2, text="选文件…", command=self._pick_ligand).pack(side="left")
        ttk.Button(r2, text="选文件夹…", command=self._pick_ligand_dir).pack(side="left", padx=(4, 0))
        ttk.Label(r2, text="(可多选；.pdbqt / .mol2 / .sdf)").pack(side="left", padx=4)

        r3 = ttk.Frame(frm)
        r3.pack(fill="x", padx=8, pady=(4, 8))
        ttk.Label(r3, text="作业名").pack(side="left")
        ttk.Entry(r3, textvariable=self.var_name, width=16).pack(side="left", padx=(4, 14))
        ttk.Label(r3, text="净电荷").pack(side="left")
        ttk.Entry(r3, textvariable=self.var_charge, width=5).pack(side="left", padx=(4, 4))
        ttk.Checkbutton(r3, text="自动判定", variable=self.var_auto_charge).pack(side="left", padx=(0, 14))
        ttk.Label(r3, text="生产时长").pack(side="left")
        ttk.Entry(r3, textvariable=self.var_ns, width=7).pack(side="left", padx=(4, 2))
        ttk.Label(r3, text="ns").pack(side="left", padx=(0, 14))
        ttk.Label(r3, text="独立重复").pack(side="left")
        ttk.Entry(r3, textvariable=self.var_reps, width=4).pack(side="left", padx=(4, 2))
        ttk.Label(r3, text="条").pack(side="left")
        ttk.Button(r3, text="导入配对表…", command=self.on_import).pack(side="left", padx=(16, 0))

        ctl = ttk.LabelFrame(self.root, text="运行控制")
        ctl.pack(fill="x", padx=10, pady=6)

        r4 = ttk.Frame(ctl)
        r4.pack(fill="x", padx=8, pady=6)
        ttk.Button(r4, text="添加蛋白+配体", command=self.on_add).pack(side="left")
        ttk.Button(r4, text="批量配对…", command=self.on_batch).pack(side="left", padx=(4, 0))
        ttk.Separator(r4, orient="vertical").pack(side="left", fill="y", padx=10)
        ttk.Label(r4, text="并行作业").pack(side="left")
        ttk.Spinbox(r4, from_=1, to=8, width=4, textvariable=self.var_workers).pack(side="left", padx=4)
        ttk.Label(r4, text="GPU 同时跑").pack(side="left")
        ttk.Spinbox(r4, from_=1, to=4, width=4, textvariable=self.var_gpu).pack(side="left", padx=4)
        ttk.Button(r4, text="开始运行", command=self.on_start).pack(side="left", padx=(14, 4))
        ttk.Button(r4, text="停止", command=self.on_stop).pack(side="left", padx=4)
        ttk.Button(r4, text="刷新", command=self.refresh).pack(side="left", padx=4)
        ttk.Button(r4, text="打开作业目录", command=self.on_open_dir).pack(side="left", padx=4)
        ttk.Button(r4, text="导出汇总", command=self.on_report).pack(side="left", padx=4)
        ttk.Button(r4, text="出图", command=self.on_figures).pack(side="left", padx=4)
        ttk.Button(r4, text="清空历史", command=self.on_clear).pack(side="left", padx=4)

        r5 = ttk.Frame(ctl)
        r5.pack(fill="x", padx=8, pady=(0, 8))
        self.var_prog = DoubleVar(value=0.0)
        self.pb = ttk.Progressbar(r5, orient="horizontal", mode="determinate",
                                  maximum=100.0, variable=self.var_prog, length=360)
        self.pb.pack(side="left")
        self.var_eta = StringVar(value="等待开始")
        ttk.Label(r5, textvariable=self.var_eta).pack(side="left", padx=10)

    def _build_table(self):
        frm = ttk.LabelFrame(self.root, text="作业列表")
        frm.pack(fill="both", expand=True, padx=10, pady=6)
        cols = ("name", "status", "prep", "check", "density", "speed", "ns")
        heads = {"name": "作业名", "status": "状态", "prep": "准备步骤",
                 "check": "体检", "density": "NPT 密度 kg/m³",
                 "speed": "ns/day", "ns": "生产 ns"}
        widths = {"name": 200, "status": 90, "prep": 90, "check": 90,
                  "density": 140, "speed": 90, "ns": 90}
        self.tree = ttk.Treeview(frm, columns=cols, show="headings", height=8)
        for c in cols:
            self.tree.heading(c, text=heads[c])
            self.tree.column(c, width=widths[c], anchor="center")
        self.tree.pack(fill="both", expand=True, side="left")
        sb = ttk.Scrollbar(frm, orient="vertical", command=self.tree.yview)
        sb.pack(side="right", fill="y")
        self.tree.configure(yscrollcommand=sb.set)

    def _build_log(self):
        frm = ttk.LabelFrame(self.root, text="运行日志")
        frm.pack(fill="both", expand=True, padx=10, pady=(6, 10))
        self.log = Text(frm, height=12, wrap="none")
        self.log.pack(fill="both", expand=True, side="left")
        sb = ttk.Scrollbar(frm, orient="vertical", command=self.log.yview)
        sb.pack(side="right", fill="y")
        self.log.configure(yscrollcommand=sb.set, state="disabled")

    # ---------------- 菜单：设置 / 关于 ----------------
    def _build_menu(self):
        bar = Menu(self.root)
        self.root.config(menu=bar)

        m_set = Menu(bar, tearoff=0)
        m_set.add_command(label="组件路径与参数…", command=self._show_settings)
        m_set.add_command(label="自动探测组件路径", command=self.on_detect_paths)
        m_set.add_command(label="运行环境自检", command=self.on_doctor)
        m_set.add_command(label="安装 AI 技能…", command=self._show_skill_install)
        m_set.add_separator()
        m_set.add_command(label="打开 config.json", command=self.on_open_config)
        bar.add_cascade(label="设置", menu=m_set)

        m_help = Menu(bar, tearoff=0)
        m_help.add_command(label="操作手册", command=self.on_open_manual)
        m_help.add_command(label="安装与环境清单",
                           command=lambda: self._open_doc("安装与环境清单.md"))
        m_help.add_command(label="安装说明",
                           command=lambda: self._open_doc("安装说明.md"))
        m_help.add_command(label="AI 操作说明（skill）",
                           command=lambda: self._open_doc("SKILL.md"))
        m_help.add_separator()
        m_help.add_command(label="关于", command=self._show_about)
        bar.add_cascade(label="帮助", menu=m_help)

    def on_detect_paths(self):
        """探测本机组件路径并写入 config.json。"""
        try:
            import autoconfig
            cfg = mdpipe.config()
            det = autoconfig.detect_all(cfg)
            if not det:
                messagebox.showwarning(
                    "自动探测",
                    "没有探测到任何组件。请确认 GROMACS 已解压、"
                    "WSL 已安装，或在「组件路径与参数」里手动填写。")
                return
            merged = autoconfig.merge_detected(cfg, det)
            mdpipe.save_json(mdpipe.CONFIG_PATH, merged)
            self.cfg = mdpipe.config()
            lines = [f"{k} = {v}" for k, v in det.items()]
            self._log("自动探测并写入 config.json：\n" + "\n".join(lines) + "\n")
            messagebox.showinfo("自动探测", "已写入 config.json：\n\n" + "\n".join(lines))
        except Exception as exc:
            messagebox.showerror("自动探测", f"探测失败：{exc}")

    def _show_settings(self):
        dlg = Toplevel(self.root)
        dlg.title("设置：组件路径")
        dlg.transient(self.root)
        dlg.resizable(False, False)
        frm = ttk.Frame(dlg, padding=12)
        frm.pack(fill="both", expand=True)

        ttk.Label(frm, text="这些值保存在 pipeline/config.json，保存后对下一次运行立即生效。",
                  foreground="#555").grid(row=0, column=0, columnspan=4, sticky="w",
                                          pady=(0, 10))
        cfg = mdpipe.config()
        vars_: dict[str, StringVar] = {}
        for i, (key, label, hint, browse) in enumerate(PATH_FIELDS, start=1):
            ttk.Label(frm, text=label).grid(row=i, column=0, sticky="e",
                                            padx=(0, 8), pady=3)
            v = StringVar(value=str(cfg.get(key, "")))
            vars_[key] = v
            ttk.Entry(frm, textvariable=v, width=52).grid(row=i, column=1,
                                                          sticky="we", pady=3)
            if browse:
                def pick(v=v, key=key):
                    if key == "gmx":
                        p = filedialog.askopenfilename(
                            title="选择 gmx.exe", parent=dlg,
                            filetypes=[("gmx.exe", "gmx.exe"), ("可执行文件", "*.exe")])
                    else:
                        p = filedialog.askdirectory(title="选择 LigPlot+ 根目录",
                                                    parent=dlg)
                    if p:
                        v.set(p)
                ttk.Button(frm, text="浏览…", width=8, command=pick).grid(
                    row=i, column=2, padx=(6, 0))
            else:
                ttk.Label(frm, text="").grid(row=i, column=2)
            ttk.Label(frm, text=hint, foreground="#777").grid(
                row=i, column=3, sticky="w", padx=(10, 0))

        ttk.Label(frm, text="路径里的正斜杠、反斜杠都可以；WSL 内路径以 /mnt/ 开头。",
                  foreground="#555").grid(row=len(PATH_FIELDS) + 1, column=0,
                                          columnspan=4, sticky="w", pady=(8, 4))

        btns = ttk.Frame(frm)
        btns.grid(row=len(PATH_FIELDS) + 2, column=0, columnspan=4, sticky="we",
                  pady=(10, 0))

        def restore():
            defaults = mdpipe.DEFAULT_CONFIG
            for key, _label, _hint, _b in PATH_FIELDS:
                vars_[key].set(str(defaults.get(key, "")))

        def save():
            new = {k: v.get().strip() for k, v in vars_.items()}
            g = new.get("gmx", "")
            if g and ("\\" in g or "/" in g) and not Path(g).exists():
                if not messagebox.askyesno(
                        "路径不存在",
                        f"{g}\n\n这个文件现在不存在，仍然保存吗？", parent=dlg):
                    return
            cfg2 = mdpipe.config()
            cfg2.update(new)
            mdpipe.save_json(mdpipe.CONFIG_PATH, cfg2)
            self.cfg = mdpipe.config()
            self._log(f"设置已保存：{mdpipe.CONFIG_PATH}\n")
            dlg.destroy()

        def detect_into_fields():
            try:
                import autoconfig
                det = autoconfig.detect_all(mdpipe.config())
            except Exception as exc:
                messagebox.showerror("自动探测", f"探测失败：{exc}", parent=dlg)
                return
            filled = 0
            for key, value in det.items():
                if key in vars_ and value:
                    vars_[key].set(value)
                    filled += 1
            if not filled:
                messagebox.showinfo("自动探测", "没有探测到可用的组件路径。", parent=dlg)

        ttk.Button(btns, text="自动探测", command=detect_into_fields).pack(side="left")
        ttk.Button(btns, text="恢复默认", command=restore).pack(side="left")
        ttk.Button(btns, text="运行环境自检", command=self.on_doctor).pack(side="left", padx=8)
        ttk.Button(btns, text="取消", command=dlg.destroy).pack(side="right")
        ttk.Button(btns, text="保存", command=save).pack(side="right", padx=(0, 6))
        dlg.grab_set()

    def _show_skill_install(self):
        """把 md-pipeline 技能装到用户选中的 AI Harness；只复制文本文件。"""
        try:
            import install_skill
        except Exception as e:
            messagebox.showerror("安装 AI 技能", f"找不到 install_skill.py：{e}")
            return

        dlg = Toplevel(self.root)
        dlg.title("安装 AI 技能（md-pipeline）")
        dlg.transient(self.root)
        dlg.geometry("980x540")

        ttk.Label(
            dlg, padding=(12, 10, 12, 4),
            text="勾选要安装的目标。检测=有 表示本机已装该 Harness；"
                 "规则文件（AGENTS.md / Trae / Qoder / Cursor）只写一段指引。"
        ).pack(anchor="w")

        body = ttk.Frame(dlg, padding=(12, 0))
        body.pack(fill="both", expand=True)
        lb = Listbox(body, selectmode="extended", width=120, height=18,
                     font=("Consolas", 9), activestyle="none")
        sb = ttk.Scrollbar(body, command=lb.yview)
        lb.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        lb.pack(fill="both", expand=True)

        reg = install_skill.registry()
        keys: list[str] = []

        def fill():
            lb.delete(0, "end")
            keys.clear()
            for k, t in reg.items():
                keys.append(k)
                has = "有" if install_skill.detected(t) else "-"
                lb.insert("end", f"{k:<11} {install_skill.status(t):<8} "
                                 f"{'检测:' + has:<7} {install_skill.dest_of(t)}")
                if install_skill.detected(t) and install_skill.status(t) == "未安装":
                    lb.selection_set("end")
            if lb.curselection():
                lb.activate(lb.curselection()[0])   # 免得首行被误画成选中
        fill()

        msg = ttk.Label(dlg, text="", padding=(12, 4), foreground="#3E8E52")
        msg.pack(anchor="w")

        btns = ttk.Frame(dlg, padding=(12, 6, 12, 12))
        btns.pack(fill="x")

        def choose_detected():
            lb.selection_clear(0, "end")
            for i, k in enumerate(keys):
                if install_skill.detected(reg[k]) and install_skill.status(reg[k]) == "未安装":
                    lb.selection_set(i)

        def do_install():
            sel = [keys[i] for i in lb.curselection()]
            if not sel:
                messagebox.showinfo("安装 AI 技能", "先选一个目标。", parent=dlg)
                return
            text = install_skill.SOURCE.read_text(encoding="utf-8")
            done, exists = [], []
            for k in sel:
                state = install_skill.install(k, text, force=False, reg=reg)
                if state == "exists":
                    exists.append(k)
                elif state in ("installed", "unchanged"):
                    done.append(k)
            if exists:
                if messagebox.askyesno(
                        "已存在不同版本",
                        "以下目标已有不同版本的技能，覆盖吗？\n\n"
                        + "\n".join(exists), parent=dlg):
                    for k in exists:
                        install_skill.install(k, text, force=True, reg=reg)
                    done += exists
            fill()
            msg.configure(text=f"完成：{'、'.join(done) if done else '无变化'}")

        ttk.Button(btns, text="选中检测到的", command=choose_detected).pack(side="left")
        ttk.Button(btns, text="刷新", command=fill).pack(side="left", padx=8)
        ttk.Button(btns, text="关闭", command=dlg.destroy).pack(side="right")
        ttk.Button(btns, text="安装所选", command=do_install).pack(side="right", padx=(0, 6))
        dlg.grab_set()

    def _show_about(self):
        dlg = Toplevel(self.root)
        dlg.title("关于")
        dlg.transient(self.root)
        dlg.geometry("640x680")
        frm = ttk.Frame(dlg, padding=10)
        frm.pack(fill="both", expand=True)
        txt = Text(frm, wrap="word", relief="flat", background="#fafafa")
        txt.pack(fill="both", expand=True)
        txt.insert("1.0", ABOUT_TEXT.format(version=APP_VERSION))
        txt.configure(state="disabled")

        bar = ttk.Frame(frm)
        bar.pack(fill="x", pady=(8, 0))
        ttk.Button(bar, text="运行环境自检", command=self.on_doctor).pack(side="left")
        ttk.Button(bar, text="关闭", command=dlg.destroy).pack(side="right")
        dlg.grab_set()

    def on_doctor(self):
        self._log("运行环境自检…\n")

        def body():
            try:
                r = subprocess.run(
                    [sys.executable, str(mdpipe.ROOT / "scripts" / "doctor.py")],
                    cwd=str(mdpipe.ROOT), capture_output=True, text=True,
                    encoding="utf-8", errors="replace", timeout=180)
                out = (r.stdout or "") + (r.stderr or "")
            except Exception as exc:  # noqa: BLE001
                out = f"自检失败：{exc}\n"
            self.root.after(0, lambda: self._log(out + "\n"))

        threading.Thread(target=body, daemon=True).start()

    def on_open_config(self):
        mdpipe.CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
        if not mdpipe.CONFIG_PATH.exists():
            mdpipe.save_json(mdpipe.CONFIG_PATH, mdpipe.DEFAULT_CONFIG)
        os.startfile(str(mdpipe.CONFIG_PATH))

    def on_open_manual(self):
        self._open_doc("分子动力学操作手册.md")

    def _open_doc(self, name: str):
        p = _find_doc(name)
        if p is None:
            messagebox.showinfo("打开文档", f"没找到 {name}\n安装目录：{mdpipe.ROOT}")
            return
        try:
            win = mdview.show(self.root, p)
            icon = _icon_path()
            if icon:
                try:
                    win.iconbitmap(default=str(icon))
                except Exception:
                    pass
        except Exception:
            # 渲染器本身出问题时的最后兜底：纯文本
            self._show_text_file(p)

    def _show_text_file(self, p: Path):
        dlg = Toplevel(self.root)
        dlg.title(p.name)
        dlg.geometry("900x700")
        txt = Text(dlg, wrap="word", font=("Consolas", 10))
        sb = ttk.Scrollbar(dlg, command=txt.yview)
        txt.configure(yscrollcommand=sb.set)
        txt.insert("1.0", p.read_text(encoding="utf-8", errors="replace"))
        txt.configure(state="disabled")
        sb.pack(side="right", fill="y")
        txt.pack(fill="both", expand=True)
        ttk.Label(dlg, text=f"文件位置：{p}", foreground="#6B7785").pack(
            anchor="w", padx=8, pady=4)

    # ---------------- 动作 ----------------
    def _pick_protein(self):
        f = filedialog.askopenfilename(title="选择蛋白结构",
                                       filetypes=[("PDB", "*.pdb *.ent"), ("全部", "*.*")])
        if f:
            self.var_protein.set(f)
            if not self.var_name.get():
                self.var_name.set(Path(f).stem)

    def _pick_complex(self):
        fs = filedialog.askopenfilenames(
            title="选择复合物 PDB（可按住 Ctrl 多选）",
            filetypes=[("PDB", "*.pdb *.ent"), ("全部", "*.*")])
        if fs:
            self.var_complex.set(" ; ".join(fs))

    def _complex_files(self):
        raw = self.var_complex.get().replace("|", ";").replace("\n", ";")
        return [c.strip().strip('"') for c in raw.split(";") if c.strip().strip('"')]

    def on_add_complex(self):
        """复合物入口：一个文件一个作业，自动拆成蛋白 + 配体。"""
        files = self._complex_files()
        if not files:
            messagebox.showwarning("缺少文件", "先用「选文件…」选一个或多个复合物 PDB")
            return
        nums = self._num_fields()
        if nums is None:
            return
        ns, reps, charge = nums

        spec = type("A", (), {})()
        spec.complex = files
        spec.name = self.var_name.get().strip() or None
        spec.lig_resname = None
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            plan, extras = mdpipe._plan_from_complex(spec)
        self._log(buf.getvalue())
        if not plan:
            messagebox.showwarning(
                "没有可登记的作业",
                "这些文件里没找到配体和蛋白共存——确认是「蛋白+配体」的复合物 PDB，"
                "而不是单独的蛋白或单独的配体。")
            return
        self._register(plan, ns, reps, charge, extras)

    def _pick_ligand(self):
        fs = filedialog.askopenfilenames(title="选择配体位姿（可按住 Ctrl 多选）",
                                         filetypes=[("对接位姿", "*.pdbqt"),
                                                    ("分子文件", "*.mol2 *.sdf"),
                                                    ("全部", "*.*")])
        if fs:
            self._set_ligands(list(fs))

    def _pick_ligand_dir(self):
        d = filedialog.askdirectory(title="选择配体文件夹")
        if d:
            self._set_ligands([str(p) for p in mdpipe.find_files([d])
                               if p.suffix.lower() in mdpipe.LIGANT_EXTS])

    def _set_ligands(self, paths: list[str]):
        self.ligand_list = paths
        if not paths:
            self.var_ligand.set("")
            return
        if len(paths) == 1:
            self.var_ligand.set(paths[0])
        else:
            self.var_ligand.set(f"{Path(paths[0]).name} 等 {len(paths)} 个文件")

    def _ligands(self) -> list[str]:
        if self.ligand_list:
            return self.ligand_list
        text = self.var_ligand.get().strip()
        if not text:
            return []
        parts = [t.strip().strip('"') for t in text.split(";") if t.strip()]
        return parts or [text]

    def _num_fields(self):
        try:
            ns = float(self.var_ns.get())
            reps = int(self.var_reps.get())
        except ValueError:
            messagebox.showwarning("参数有误", "生产时长和重复数必须是数字")
            return None
        charge = None
        if not self.var_auto_charge.get():
            try:
                charge = int(self.var_charge.get())
            except ValueError:
                messagebox.showwarning("参数有误", "净电荷必须是整数")
                return None
        return ns, reps, charge

    def _register(self, plan, ns, reps, charge, extras=None):
        """按 (name, protein, ligand) 清单逐个登记；charge=None 时逐个自动判定。"""
        extras = extras or {}
        added = 0
        for name, protein, ligand in plan:
            ch = charge
            if ch is None:
                ch = mdpipe.detect_charge(self.cfg, protein, ligand)
                if ch is None:
                    ch = 0
                    self._log(f"{name}: 未能自动判定净电荷，按 0 处理——请在 job.json 里核对\n")
                else:
                    self._log(f"{name}: 净电荷 {ch}\n")
            meta = {"name": name, "protein": str(Path(protein).resolve()),
                    "ligand": str(Path(ligand).resolve()), "charge": int(ch),
                    "resname": "LIG", "model": 1, "replicates": reps, "prod_ns": ns}
            meta.update(extras.get(name, {}))
            if mdpipe.add_job(meta):
                added += 1
        self._log(f"本次登记 {added} 个作业\n")
        self.refresh()

    def on_add(self):
        protein = self.var_protein.get().strip()
        ligands = self._ligands()

        # 只选了蛋白：如果它其实是复合物，拆开继续
        if protein and not ligands:
            lig_res = mdpipe.detect_complex(protein)
            if lig_res:
                if not messagebox.askyesno(
                        "检测到复合物",
                        f"{Path(protein).name} 里有配体残基 {lig_res}，是复合物。\n\n"
                        f"这类文件请用上面专门的「复合物 PDB」那一行。\n"
                        f"现在就按复合物拆开继续添加吗？"):
                    return
                name = mdpipe.sanitize_name(self.var_name.get().strip() or Path(protein).stem)
                try:
                    prot, ligf = mdpipe.split_complex_job(protein, name, None)
                except Exception as e:
                    messagebox.showerror("拆分失败", str(e))
                    return
                self._log(f"{Path(protein).name}: 已拆出配体 {lig_res}\n")
                nums = self._num_fields()
                if nums is None:
                    return
                ns, reps, charge = nums
                self._register([(name, prot, ligf)], ns, reps, charge)
                return

        if not protein or not ligands:
            messagebox.showwarning("缺少文件",
                                   "请在「蛋白 PDB」选蛋白、在「配体位姿」选配体；\n"
                                   "如果手上是蛋白+配体合在一起的复合物，用「复合物 PDB」那一行。")
            return
        nums = self._num_fields()
        if nums is None:
            return
        ns, reps, charge = nums

        proteins = [Path(protein)]
        ligs = [Path(l) for l in ligands]
        plan = mdpipe.plan_jobs(proteins, ligs, self.var_name.get().strip() or None)
        self._register(plan, ns, reps, charge)

    def on_batch(self):
        """多蛋白 × 多配体：选两个文件夹，做交叉组合。"""
        pdir = filedialog.askdirectory(title="第 1 步：选择蛋白文件夹")
        if not pdir:
            return
        ldir = filedialog.askdirectory(title="第 2 步：选择配体文件夹")
        if not ldir:
            return
        proteins = [p for p in mdpipe.find_files([pdir])
                    if p.suffix.lower() in mdpipe.PROTEIN_EXTS]
        ligands = [p for p in mdpipe.find_files([ldir])
                   if p.suffix.lower() in mdpipe.LIGANT_EXTS]
        if not proteins or not ligands:
            messagebox.showwarning("没有文件", "选中的文件夹里没有可用的蛋白或配体文件")
            return
        all_pairs = [(p, l) for p in proteins for l in ligands]
        picked = PairDialog(self.root, all_pairs,
                            title=f"选择配对（{len(proteins)} 蛋白 × {len(ligands)} 配体）").show()
        if not picked:
            self._log("没有选择任何配对\n")
            return
        plan = mdpipe.plan_pairs(picked)
        nums = self._num_fields()
        if nums is None:
            return
        ns, reps, charge = nums
        self._register(plan, ns, reps, charge)

    def on_import(self):
        """从 CSV 配对表登记：一行一个蛋白–配体对。"""
        f = filedialog.askopenfilename(
            title="选择配对表（需要 protein 和 ligand 两列）",
            filetypes=[("表格", "*.csv *.tsv *.txt"), ("全部", "*.*")])
        if not f:
            return
        table = Path(f).resolve()
        rows = mdpipe.resolve_table_rows(mdpipe.load_pairs_table(table), table.parent)
        if not rows:
            messagebox.showwarning("表格为空", "表里至少要有一行 protein 和 ligand")
            return
        nums = self._num_fields()
        if nums is None:
            return
        ns, reps, charge = nums
        added = 0
        for row in rows:
            ch = row.get("charge")
            ch = int(ch) if ch not in ("", None) else charge
            if ch is None:
                ch = mdpipe.detect_charge(self.cfg, Path(row["protein"]), Path(row["ligand"]))
                if ch is None:
                    ch = 0
                    self._log(f"{Path(row['ligand']).name}: 净电荷按 0 处理\n")
            name = row.get("name") or mdpipe.sanitize_name(
                f"{Path(row['protein']).stem}__{Path(row['ligand']).stem}")
            meta = {
                "name": name, "protein": row["protein"], "ligand": row["ligand"],
                "charge": int(ch),
                "resname": row.get("resname") or "LIG",
                "model": int(row.get("model") or 1),
                "replicates": int(row.get("replicates") or reps),
                "prod_ns": float(row.get("ns") or ns),
            }
            if mdpipe.add_job(meta):
                added += 1
        self._log(f"从配对表登记 {added} 个作业\n")
        self.refresh()

    def on_start(self):
        if self.running:
            messagebox.showinfo("正在运行", "队列已经在跑了")
            return
        q = mdpipe.queue()
        todo = [n for n in q["jobs"] if mdpipe.job_meta(n).get("status") != "done"]
        if not todo:
            messagebox.showinfo("没有待跑作业", "队列为空，或作业都已完成")
            return

        self.running = True
        mdpipe.STOP.clear()
        workers = int(self.var_workers.get() or 2)
        gpu = int(self.var_gpu.get() or 2)
        self._log(f"开始：{len(todo)} 个作业，并行 {workers}，GPU 同时跑 {gpu}\n")

        def body():
            old = sys.stdout
            sys.stdout = QueueWriter(self.logq)
            try:
                sem = threading.Semaphore(gpu)
                with ThreadPoolExecutor(max_workers=workers) as pool:
                    futs = [pool.submit(mdpipe.worker, self.cfg, n, sem, None) for n in todo]
                    for f in as_completed(futs):
                        f.result()
            except Exception:
                self.logq.put(traceback.format_exc())
            finally:
                sys.stdout = old
                self.running = False
                self.logq.put("__DONE__")

        self.worker_thread = threading.Thread(target=body, daemon=True)
        self.worker_thread.start()
        self.root.after(500, self._tick)

    def on_stop(self):
        if not self.running:
            return
        mdpipe.STOP.set()
        self._log("已请求停止：当前这一步结束后停下，已完成的步骤会保留\n")

    def on_open_dir(self):
        sel = self.tree.selection()
        p = mdpipe.job_dir(sel[0]) if sel else mdpipe.PIPE
        if p.exists():
            os.startfile(str(p))

    def on_report(self):
        mdpipe.cmd_report(type("A", (), {})())
        self._log(f"汇总已写出：{mdpipe.PIPE / 'summary.csv'}\n")

    def on_figures(self):
        """选一个作业出六联图；选多个出体系对比图。"""
        names = [self.tree.item(i, "values")[0] for i in self.tree.selection()]
        names = [n for n in names if mdpipe.job_dir(n).exists()]
        if not names:
            messagebox.showinfo("出图", "先在下面的列表里选中作业（按住 Ctrl 可多选）")
            return
        self._log(f"出图：{'、'.join(names)}\n")

        def body():
            try:
                import make_figures
                if len(names) == 1:
                    out = make_figures.plot_job(mdpipe.job_dir(names[0]), names[0])
                else:
                    out = make_figures.plot_compare(
                        [mdpipe.job_dir(n) for n in names], names,
                        mdpipe.PIPE / "figures")
                msg = f"图已生成：{out}\n" if out else "没有可画的数据，先确认分析步骤跑完\n"
            except Exception as exc:  # noqa: BLE001
                msg = f"出图失败：{exc}\n"
            self.root.after(0, lambda: self._log(msg))

        threading.Thread(target=body, daemon=True).start()

    def on_clear(self):
        """把已完成的作业移出队列；问清楚要不要连目录一起删。"""
        q = mdpipe.queue()
        done_jobs = [n for n in q["jobs"]
                     if mdpipe.job_meta(n).get("status") == "done"]
        if not done_jobs:
            messagebox.showinfo("清空历史", "没有已完成的作业")
            return
        ans = messagebox.askyesnocancel(
            "清空历史",
            f"有 {len(done_jobs)} 个已完成的作业要移出队列。\n\n"
            "「是」= 移出队列并删除作业目录（轨迹、分析和图一起删）\n"
            "「否」= 只移出队列，目录留在硬盘上\n"
            "「取消」= 什么都不做")
        if ans is None:
            return
        q["jobs"] = [n for n in q["jobs"] if n not in done_jobs]
        mdpipe.save_queue(q)
        if ans:
            for n in done_jobs:
                shutil.rmtree(mdpipe.job_dir(n), ignore_errors=True)
        self._log("清空历史：移出 %d 个作业%s\n"
                  % (len(done_jobs), "，目录已删除" if ans else "，目录保留"))
        self.refresh()

    # ---------------- 刷新 ----------------
    def _tick(self):
        self.refresh()
        if self.running:
            self.root.after(2000, self._tick)

    def refresh(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        q = mdpipe.queue()
        for name in q["jobs"]:
            meta = mdpipe.job_meta(name)
            if not meta:
                continue
            prep = sum(1 for st in mdpipe.PREP_STEPS if mdpipe.done(name, st))
            m = meta.get("metrics", {}) or {}
            dens = m.get("npt_density")
            spd = m.get("ns_per_day")
            chk = meta.get("structure_check") or {}
            label = chk.get("verdict", "-")
            if chk.get("repaired"):
                label += "·已补"
            self.tree.insert("", "end", iid=name, values=(
                name, meta.get("status", "?"), f"{prep}/8", label,
                f"{dens:.1f}" if isinstance(dens, (int, float)) else "-",
                f"{spd:.0f}" if isinstance(spd, (int, float)) else "-",
                meta.get("prod_ns", "-")))
        self._update_progress()

    def _update_progress(self):
        """进度条按"已完成步骤/总步骤"算；剩余时间按剩余 ns ÷ 实测吞吐估。"""
        try:
            p = mdpipe.progress()
        except Exception:  # noqa: BLE001
            return
        self.var_prog.set(p["pct"])
        if p["jobs_total"] == 0:
            self.var_eta.set("队列为空")
            return
        if p["jobs_done"] == p["jobs_total"]:
            self.var_prog.set(100.0)
            self.var_eta.set(f"全部完成（{p['jobs_total']} 个作业）")
            return
        h = p["eta_hours"]
        eta = f"{h:.1f} 小时" if h < 48 else f"{h / 24:.1f} 天"
        self.var_eta.set(
            f"完成 {p['jobs_done']}/{p['jobs_total']} 个作业 · "
            f"步骤 {p['done_steps']}/{p['total_steps']} · "
            f"剩余 {p['remaining_ns']:.0f} ns · 约 {eta}"
            f"（按 {p['speed']:.0f} ns/day 估）")

    def _drain_log(self):
        try:
            while True:
                text = self.logq.get_nowait()
                if text == "__DONE__":
                    self.running = False
                    self.refresh()
                    continue
                self._log(text)
        except queue.Empty:
            pass
        self.root.after(200, self._drain_log)

    def _log(self, text):
        self.log.configure(state="normal")
        self.log.insert("end", text)
        self.log.see("end")
        self.log.configure(state="disabled")


def main():
    try:
        root = Tk()
        App(root)
        root.mainloop()
    except Exception:
        # 用 pythonw 从快捷方式启动时没有控制台，出错要弹窗并留日志
        text = traceback.format_exc()
        log = mdpipe.PIPE / "gui_error.log"
        try:
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(text, encoding="utf-8")
        except Exception:
            pass
        try:
            messagebox.showerror("界面启动失败", f"{text}\n\n详情已写入 {log}")
        except Exception:
            pass


if __name__ == "__main__":
    main()
