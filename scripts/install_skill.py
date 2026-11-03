#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""install_skill.py —— 把 md-pipeline 技能装到各种 AI Harness 认识的位置。

    python scripts/install_skill.py --list
    python scripts/install_skill.py --target clawsgo dsh workbuddy
    python scripts/install_skill.py --target detected        # 只装到本机检测到的
    python scripts/install_skill.py --target agents          # 往项目 AGENTS.md 追加
    python scripts/install_skill.py --print

两种落点：

* skill 目录 —— `<根目录>/md-pipeline/SKILL.md`，主流 Harness 都是这个约定
  （ClawsGO、DeepSeek Harness、WorkBuddy、Modao、Codex、Claude Code、~/.agents）。
* 规则文件 —— 单个 Markdown，写一段精简规则并指向完整 SKILL.md
  （AGENTS.md、Trae、Qoder、Cursor、Copilot）。

只是复制文本，不执行、不联网、不改注册表。
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "skills" / "md-pipeline" / "SKILL.md"

AGENTS_BEGIN = "<!-- md-pipeline:begin -->"
AGENTS_END = "<!-- md-pipeline:end -->"


@dataclass
class Target:
    key: str
    label: str
    kind: str          # skill | rule
    path: Path         # skill: 根目录；rule: 文件
    detect: Path       # 用来判断“这个 Harness 装没装”的目录/文件
    note: str = ""


def registry(home: Path | None = None, cwd: Path | None = None) -> dict[str, Target]:
    home = home or Path.home()
    cwd = cwd or Path.cwd()
    rows = [
        # —— skill 目录：<根>/md-pipeline/SKILL.md
        ("clawsgo", "ClawsGO 智能体", "skill",
         home / ".clawsgo" / "skills", home / ".clawsgo", "安装包默认目标"),
        ("dsh", "DeepSeek Harness", "skill",
         home / ".dsh" / "skills", home / ".dsh", "本机技能库最大"),
        ("workbuddy", "WorkBuddy", "skill",
         home / ".workbuddy" / "skills", home / ".workbuddy", ""),
        ("modao", "Modao（pi）", "skill",
         home / ".modao" / "agents" / "pi" / "skills", home / ".modao", ""),
        ("agentshome", "通用 ~/.agents", "skill",
         home / ".agents" / "skills", home / ".agents", "多个工具共用"),
        ("codex", "Codex CLI", "skill",
         home / ".codex" / "skills", home / ".codex", ""),
        ("claude", "Claude Code", "skill",
         home / ".claude" / "skills", home / ".claude", ""),
        # —— 项目级 skill 目录
        ("project", "本项目（ClawsGO）", "skill",
         cwd / ".clawsgo" / "skills", cwd / ".clawsgo", "只影响当前项目"),
        ("pclaude", "本项目（Claude）", "skill",
         cwd / ".claude" / "skills", cwd / ".claude", "只影响当前项目"),
        # —— 规则文件：写一个精简规则并指向 SKILL.md
        ("agents", "AGENTS.md", "rule",
         cwd / "AGENTS.md", cwd / "AGENTS.md", "通用智能体入口"),
        ("trae", "Trae CN（项目规则）", "rule",
         cwd / ".trae" / "rules" / "md-pipeline.md", cwd / ".trae", ""),
        ("qoder", "Qoder（项目规则）", "rule",
         cwd / ".qoder" / "rules" / "md-pipeline.md", cwd / ".qoder", ""),
        ("cursor", "Cursor（项目规则）", "rule",
         cwd / ".cursor" / "rules" / "md-pipeline.mdc", cwd / ".cursor", ""),
        ("copilot", "Copilot（项目指令）", "rule",
         cwd / ".github" / "instructions" / "md-pipeline.instructions.md",
         cwd / ".github", ""),
    ]
    return {r[0]: Target(*r) for r in rows}


def rule_body(source: Path = SOURCE) -> str:
    return f"""{AGENTS_BEGIN}
## 分子动力学（mdpipe）

本机装有 ClawsGO MD 流水线。要跑 MD、看轨迹、算 MM-PBSA 时，先读
`{source}` 的完整说明，再动手。

顺序：`mdpipe doctor` → `mdpipe status` → `mdpipe run --job …` →
`mdpipe analysis` / `mdpipe mmpbsa` / `mdpipe interactions`。

硬约束：不要在有作业 running 时再起 mdrun；分析轨迹用 `md_center.xtc` /
`md_prot.xtc`，禁止 `-pbc nojump`；`pipeline/jobs` 是数据，清理用 `mdpipe archive`。
{AGENTS_END}
"""


def install_one(base: Path, text: str, force: bool = False) -> str:
    """skill 目录安装；返回 installed / unchanged / exists。"""
    dest = base / "md-pipeline" / "SKILL.md"
    if dest.is_file():
        if dest.read_text(encoding="utf-8", errors="replace") == text:
            return "unchanged"
        if not force:
            return "exists"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(text, encoding="utf-8")
    return "installed"


def install_rule(path: Path, text: str, force: bool = False,
                 cursor: bool = False) -> str:
    """规则文件安装；带标记的段落可重复写入。"""
    old = path.read_text(encoding="utf-8") if path.is_file() else ""
    if AGENTS_BEGIN in old and AGENTS_END in old:
        # 已装过：只替换标记段，不动文件其余内容和 Cursor 的 frontmatter
        new = old.replace(
            old[old.index(AGENTS_BEGIN):old.index(AGENTS_END) + len(AGENTS_END)],
            text.strip())
        if new == old:
            return "unchanged"
        if not force:
            return "exists"
    else:
        if old and not force:
            return "exists"
        new = ((old.rstrip() + "\n\n") if old.strip() else "") + text
        if cursor:
            new = ("---\ndescription: GROMACS 分子动力学流水线 mdpipe 的操作说明\n"
                   "alwaysApply: false\n---\n\n" + new)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(new, encoding="utf-8")
    return "installed"


def install_agents(path: Path, force: bool = False) -> str:
    """兼容旧接口：等价于往 AGENTS.md 写规则段。"""
    return install_rule(path, rule_body(), force)


def status(t: Target) -> str:
    if t.kind == "skill":
        return "已安装" if (t.path / "md-pipeline" / "SKILL.md").is_file() else "未安装"
    if t.path.is_file():
        return "已写入" if AGENTS_BEGIN in t.path.read_text(
            encoding="utf-8", errors="replace") else "文件存在（无本技能）"
    return "未写入"


def detected(t: Target) -> bool:
    return t.detect.exists()


def dest_of(t: Target) -> Path:
    return t.path if t.kind == "rule" else t.path / "md-pipeline" / "SKILL.md"


def install(key: str, text: str, force: bool = False,
            reg: dict[str, Target] | None = None) -> str:
    reg = reg or registry()
    t = reg[key]
    if t.kind == "skill":
        return install_one(t.path, text, force)
    return install_rule(t.path, text, force, cursor=(t.key == "cursor"))


def main(argv: list[str] | None = None) -> int:
    reg = registry()
    ap = argparse.ArgumentParser(description="安装 md-pipeline 技能到 AI Harness")
    ap.add_argument("--target", nargs="+", default=[],
                    choices=[*reg.keys(), "detected"],
                    help="安装目标；detected=只装到本机检测到的 Harness")
    ap.add_argument("--list", action="store_true", help="列出各目标的状态")
    ap.add_argument("--print", dest="dump", action="store_true", help="打印技能内容")
    ap.add_argument("--force", action="store_true", help="覆盖已存在的不同版本")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    if not SOURCE.is_file():
        print(f"找不到技能文件：{SOURCE}", file=sys.stderr)
        return 2
    text = SOURCE.read_text(encoding="utf-8")

    if a.dump:
        print(text)
        return 0

    if a.list or not a.target:
        print(f"技能源文件：{SOURCE}\n")
        print(f"{'目标':<11}{'状态':<8}{'检测':<6}路径")
        for k, t in reg.items():
            mark = "有" if detected(t) else "-"
            print(f"{k:<11}{status(t):<8}{mark:<6}{dest_of(t)}")
        print("\n用 --target <名字...> 安装；detected 表示只装到“检测=有”的条目。")
        return 0

    names = list(reg.keys()) if "detected" in a.target else a.target
    if "detected" in a.target:
        names = [k for k, t in reg.items() if detected(t)]

    rc = 0
    for k in names:
        state = install(k, text, a.force, reg)
        if not a.quiet:
            note = {"installed": "已安装", "unchanged": "已是最新",
                    "exists": "已存在不同版本（加 --force 覆盖）"}[state]
            print(f"{k:<11}{note:<22}{dest_of(reg[k])}")
        if state == "exists":
            rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
