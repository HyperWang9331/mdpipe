#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""批量准备受体：体检 → 尝试补链 → 复查 → 出汇总表。

一次处理一整个文件夹的受体 PDB。每个文件的处理链路是：

    check_protein  体检
      ├─ PASS ──────────────→ 直接复制到输出目录
      ├─ WARN ──────────────→ complete_protein 补侧链重原子 → 复查
      └─ FAIL ──┬─ 可补 ───→ complete_protein 补 loop → 复查
                └─ 补不动 ─→ 原样保留，汇总表里标明原因

补不动的典型情况：PDBFixer 只认 SEQRES 建缺口，而 PyMOL / AutoDock 处理过的
`*_clean.pdb` 把 SEQRES 删光了。这时用 `--seqres-dir` 指向原版 PDB 所在目录，
脚本会按文件名前缀自动借（AKT1_clean.pdb → AKT1.pdb）。

资源：每个文件单线程、BelowNormal、不碰 GPU，可以和正在跑的 MD 并行。

用法：
    python scripts/prepare_receptors.py 受体目录/ --out 准备好/
    python scripts/prepare_receptors.py 受体目录/ --out 准备好/ --seqres-dir 原始PDB目录/
    python scripts/prepare_receptors.py a.pdb b.pdb --out 准备好/ --report 汇总.md
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_protein import analyse  # noqa: E402


def _say(*a):
    print(*a, flush=True)


def find_mdprep_python() -> Path | None:
    """找装 openmm/pdbfixer 的那个解释器。"""
    cands = [Path.home() / "miniconda3" / "envs" / "mdprep" / "python.exe",
             Path.home() / "anaconda3" / "envs" / "mdprep" / "python.exe"]
    for p in cands:
        if p.exists():
            try:
                r = subprocess.run([str(p), "-c", "import openmm, pdbfixer"],
                                   capture_output=True, timeout=60)
                if r.returncode == 0:
                    return p
            except Exception:
                pass
    return None


def find_seqres_donor(target: Path, seqres_dir: Path | None) -> Path | None:
    """给 `X_clean.pdb` 找同名的原版 `X.pdb`。"""
    if not seqres_dir or not seqres_dir.is_dir():
        return None
    stem = target.stem
    for suffix in ("_clean", "_prepared", "_fixed", "_receptor"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    for cand in (seqres_dir / f"{stem}.pdb", seqres_dir / f"{stem}.cif",
                 seqres_dir / f"{stem}.mmcif"):
        if cand.exists() and cand != target:
            try:
                if any(ln.startswith("SEQRES")
                       for ln in cand.read_text(errors="replace").splitlines()[:2000]):
                    return cand
            except Exception:
                pass
    return None


def run_fix(py: Path, src: Path, dst: Path, donor: Path | None, args) -> tuple[bool, str]:
    cmd = [str(py), str(Path(__file__).with_name("complete_protein.py")),
           str(src), "--out", str(dst),
           "--max-loop", str(args.max_loop)]
    if donor:
        cmd += ["--seqres-from", str(donor)]
    if args.build_termini:
        cmd.append("--build-termini")
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800,
                       encoding="utf-8", errors="replace")
    out = (r.stdout or "") + (r.stderr or "")
    return r.returncode == 0, out.strip().splitlines()[-1] if out.strip() else ""


def summarize(rep: dict) -> tuple[str, str]:
    errs = [i for i in rep["issues"] if i["level"] == "ERROR"]
    gaps = sum(len(c["gaps_numbering"]) for c in rep["chains"])
    if errs:
        detail = "；".join(i["what"] for i in errs[:3])
        if len(errs) > 3:
            detail += f"…（共 {len(errs)} 条）"
        return rep["verdict"], detail
    return rep["verdict"], f"{gaps} 处缺口" if gaps else "完整"


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    ap = argparse.ArgumentParser(description="批量体检 + 补链 + 复查受体结构")
    ap.add_argument("inputs", nargs="+", help="PDB 文件，或含 PDB 的目录")
    ap.add_argument("--out", required=True, help="整理好的结构放这里")
    ap.add_argument("--seqres-dir", help="原版 PDB 目录，给 *_clean.pdb 借 SEQRES")
    ap.add_argument("--report", help="汇总表写到这里（Markdown）；默认写到 --out 下")
    ap.add_argument("--csv", help="另外写一份 CSV")
    ap.add_argument("--no-fix", action="store_true", help="只体检不补")
    ap.add_argument("--max-loop", type=int, default=15,
                    help="内部缺口最多补多少残基（默认 15）")
    ap.add_argument("--build-termini", action="store_true", help="首尾悬空段也补")
    a = ap.parse_args(argv)

    files: list[Path] = []
    for tok in a.inputs:
        p = Path(tok)
        if p.is_dir():
            files += sorted([q for q in p.iterdir()
                             if q.suffix.lower() in (".pdb", ".cif", ".mmcif")])
        elif p.exists():
            files.append(p)
    if not files:
        sys.exit("没找到任何 PDB 文件")

    outdir = Path(a.out)
    outdir.mkdir(parents=True, exist_ok=True)

    py = None if a.no_fix else find_mdprep_python()
    if not a.no_fix and not py:
        _say("! 没找到装 openmm/pdbfixer 的 mdprep 环境，本次只体检不补链。")
        _say("  建环境：conda create -n mdprep -c conda-forge python=3.11 openmm pdbfixer -y")

    seqdir = Path(a.seqres_dir) if a.seqres_dir else None
    rows = []

    for i, f in enumerate(files, 1):
        _say(f"[{i}/{len(files)}] {f.name}")
        try:
            before = analyse(f)
        except Exception as e:
            _say(f"    解析失败：{e}")
            rows.append({"file": f.name, "before": "ERROR", "after": "-",
                         "detail": str(e), "output": ""})
            continue

        v0, d0 = summarize(before)
        _say(f"    体检：{v0} — {d0}")

        dst = outdir / f"{f.stem}_fixed.pdb"
        v1, d1, note = v0, d0, ""

        if v0 == "PASS" or a.no_fix or not py:
            try:
                shutil.copy2(f, dst)
                _say(f"    复制 → {dst.name}")
            except Exception as e:
                _say(f"    复制失败：{e}")
        else:
            donor = find_seqres_donor(f, seqdir)
            if donor:
                _say(f"    借 SEQRES：{donor.name}")
            ok, tail = run_fix(py, f, dst, donor, a)
            if not ok and not dst.exists():
                _say(f"    补链失败：{tail}")
                v1, d1, note = v0, d0, "补链失败，原样保留"
                shutil.copy2(f, dst)
            else:
                after = analyse(dst)
                v1, d1 = summarize(after)
                _say(f"    复查：{v1} — {d1}")
                if after["verdict"] == "FAIL":
                    note = "补完仍有缺口，需手工或 MODELLER/SWISS-MODEL"

        rows.append({"file": f.name, "before": v0, "after": v1,
                     "detail": d1, "note": note,
                     "output": dst.name if dst.exists() else ""})

    # ---- 汇总
    md = ["# 受体准备汇总", "",
          "| 文件 | 体检 | 处理后 | 说明 | 输出 |",
          "|---|---|---|---|---|"]
    for r in rows:
        md.append(f'| {r["file"]} | {r["before"]} | {r["after"]} | '
                  f'{r.get("note") or r["detail"]} | {r["output"] or "—"} |')
    n_ok = sum(1 for r in rows if r["after"] == "PASS")
    n_warn = sum(1 for r in rows if r["after"] == "WARN")
    n_fail = sum(1 for r in rows if r["after"] == "FAIL")
    md += ["", f"共 {len(rows)} 个：PASS {n_ok}，WARN {n_warn}，FAIL {n_fail}。", ""]
    if n_fail:
        md += ["FAIL 的不能直接进 MD：缺口处 pdb2gmx 会插带电末端，蛋白内部凭空多"
               "一对 NH3+/COO−。要么换更长/更完整的结构，要么用 MODELLER 或 "
               "SWISS-MODEL 建模补齐。", ""]

    rpt = Path(a.report) if a.report else outdir / "受体准备汇总.md"
    rpt.write_text("\n".join(md), encoding="utf-8")
    _say(f"\n汇总 → {rpt}")
    if a.csv:
        with Path(a.csv).open("w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=["file", "before", "after", "detail", "note", "output"])
            w.writeheader()
            for r in rows:
                w.writerow(r)
        _say(f"CSV → {a.csv}")

    for r in rows:
        _say(f'  {r["after"]:5s}  {r["file"]}')
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
