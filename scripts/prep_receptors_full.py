#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从晶体结构一路做到能进 pdb2gmx：体检 → MODELLER 补 loop → PDBFixer 补侧链 → MD 清理。

单个工具各管一段，直接拼起来会卡在几个地方，这里按顺序串好：

  1. check_protein          先看缺什么。
  2. fill_loops_modeller    只补 loop，晶体已有残基一个原子都不动。必须先做——
                            PDBFixer 建的 loop 几何差，而且补完再跑 MODELLER 就
                            看不出哪儿是原来缺的了。
  3. complete_protein       补侧链重原子（晶体里常只有 N/CA/CB 的残基），顺带
                            真空约束最小化一下。MODELLER 不补侧链。
  4. md_ready               磷酸化残基降级成母体、去掉共晶小分子和水，否则
                            pdb2gmx 报 Residue 'TPO' not found。
  5. 复检，出汇总表。

哪一步没装/用不了就跳过并记一笔，不会中断整批。

用法：
    python scripts/prep_receptors_full.py 受体目录/ --out 准备好/
    python scripts/prep_receptors_full.py a.pdb b.pdb --out 准备好/ --models 3
    python scripts/prep_receptors_full.py 受体目录/ --out 准备好/ --no-modeller
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_protein import analyse  # noqa: E402
from complete_protein import _inject_seqres  # noqa: E402
from fill_loops_modeller import find_modeller_python  # noqa: E402
from prepare_receptors import find_seqres_donor, find_mdprep_python  # noqa: E402

HERE = Path(__file__).resolve().parent
ORDER = {"FAIL": 0, "WARN": 1, "PASS": 2}


def _say(*a):
    print(*a, flush=True)


def _run(cmd: list[str], timeout: int = 3600) -> tuple[int, str]:
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                       encoding="utf-8", errors="replace")
    out = ((r.stdout or "") + (r.stderr or "")).strip()
    return r.returncode, out


def _verdict(path: Path) -> tuple[str, str]:
    rep = analyse(path)
    errs = [i["what"] for i in rep["issues"] if i["level"] == "ERROR"]
    return rep["verdict"], ("；".join(errs[:2]) if errs else "结构完整")


def _seqres_ok(path: Path) -> bool:
    if path.suffix.lower() != ".pdb":
        return False
    try:
        return any(ln.startswith("SEQRES") for ln in
                   path.read_text(errors="replace").splitlines()[:2000])
    except Exception:
        return False


def _find_fasta(target: Path, fasta_dir: Path | None) -> Path | None:
    """给 `X_clean.pdb` 找同名的 `X.fasta`。"""
    if not fasta_dir or not fasta_dir.is_dir():
        return None
    stem = target.stem
    for suffix in ("_clean", "_prepared", "_fixed", "_receptor"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    for ext in (".fasta", ".fa", ".seq"):
        cand = fasta_dir / f"{stem}{ext}"
        if cand.exists():
            return cand
    return None


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    ap = argparse.ArgumentParser(description="受体四步准备：体检→补loop→补侧链→清理")
    ap.add_argument("inputs", nargs="+", help="PDB/CIF 文件或目录")
    ap.add_argument("--out", required=True, help="成品放这里")
    ap.add_argument("--seqres-dir", help="原版 PDB 目录，给 *_clean.pdb 借 SEQRES")
    ap.add_argument("--seq-dir",
                    help="FASTA 序列目录（{原名}.fasta）。文件里连 SEQRES 都没有时用，"
                         "序列可从 UniProt 下载")
    ap.add_argument("--models", type=int, default=3, help="MODELLER 每个缺口建几个模型")
    ap.add_argument("--no-modeller", action="store_true", help="跳过 MODELLER 补 loop")
    ap.add_argument("--report", help="汇总表路径")
    ap.add_argument("--csv")
    a = ap.parse_args(argv)

    files: list[Path] = []
    for tok in a.inputs:
        p = Path(tok)
        if p.is_dir():
            files += sorted(q for q in p.iterdir()
                            if q.suffix.lower() in (".pdb", ".cif", ".mmcif"))
        elif p.exists():
            files.append(p)
    if not files:
        sys.exit("没找到任何结构文件")

    outdir = Path(a.out)
    outdir.mkdir(parents=True, exist_ok=True)
    modpy = None if a.no_modeller else find_modeller_python()
    fixpy = find_mdprep_python()
    if not modpy:
        _say("! 没有 MODELLER 环境，本次只补侧链（缺的 loop 补不了）")
    if not fixpy:
        _say("! 没有 mdprep 环境，跳过侧链补齐")

    seqdir = Path(a.seqres_dir) if a.seqres_dir else None
    fastadir = Path(a.seq_dir) if a.seq_dir else None
    rows = []

    for i, f in enumerate(files, 1):
        _say(f"[{i}/{len(files)}] {f.name}")
        try:
            v0, d0 = _verdict(f)
        except Exception as e:
            _say(f"    解析失败：{e}")
            rows.append({"file": f.name, "v0": "ERROR", "v1": "-", "v2": "-",
                         "v3": "-", "note": str(e), "output": ""})
            continue
        _say(f"    体检：{v0} — {d0}")

        tmp = Path(tempfile.mkdtemp(prefix="prep_", dir=str(outdir)))
        cur, v1, v2, v3, notes = f, "-", "-", "-", []
        dst = outdir / f"{f.stem}_ready.pdb"
        try:
            if v0 == "PASS":
                pass
            else:
                src = f
                donor = find_seqres_donor(f, seqdir)
                if donor and f.suffix.lower() == ".pdb":
                    try:
                        src = _inject_seqres(f, donor)
                        notes.append(f"借 SEQRES：{donor.name}")
                    except SystemExit as e:
                        notes.append(str(e))

                # ---- 2. MODELLER 补 loop
                fasta = _find_fasta(f, fastadir)
                if modpy and (_seqres_ok(src) or fasta):
                    tmp_out = tmp / "modeled.pdb"
                    cmd = [sys.executable, str(HERE / "fill_loops_modeller.py"),
                           str(src), "--out", str(tmp_out), "--models", str(a.models)]
                    if fasta:
                        cmd += ["--seq", str(fasta)]
                        notes.append(f"借序列：{fasta.name}")
                    rc, log = _run(cmd, timeout=7200)
                    if rc == 0 and tmp_out.exists():
                        v, d = _verdict(tmp_out)
                        v1 = v
                        if ORDER[v] >= ORDER[v0]:
                            cur = tmp_out
                        if v != "PASS":
                            notes.append(f"补 loop 后仍 {v}：{d}")
                    else:
                        tail = log.splitlines()[-1] if log else "未知错误"
                        notes.append(f"MODELLER 没成功（{tail}）")
                elif f.suffix.lower() == ".pdb":
                    notes.append("没有 SEQRES，跳过 MODELLER（加 --seq-dir 提供 FASTA 序列）")

                # ---- 3. PDBFixer 补侧链
                if fixpy:
                    tmp_out = tmp / "fixed2.pdb"
                    # 当前结构已经没有缺口时只补侧链：PDBFixer 的 SEQRES 对齐碰到
                    # 修饰残基/重排编号会误判缺口，硬塞几段进去反而把接缝拉成 12 Å
                    cur_v = _verdict(cur)[0]
                    maxloop = "0" if cur_v in ("PASS", "WARN") else "15"
                    cmd = [str(fixpy), str(HERE / "complete_protein.py"),
                           str(cur), "--out", str(tmp_out), "--max-loop", maxloop]
                    rc, log = _run(cmd, timeout=3600)
                    if rc == 0 and tmp_out.exists():
                        v, d = _verdict(tmp_out)
                        v2 = v
                        if ORDER[v] >= ORDER[_verdict(cur)[0]]:
                            cur = tmp_out
                        if v != "PASS":
                            notes.append(f"补侧链后仍 {v}：{d}")
                    else:
                        tail = log.splitlines()[-1] if log else "未知错误"
                        notes.append(f"补侧链没成功（{tail}）")

            # ---- 4. MD 清理
            rc, log = _run([sys.executable, str(HERE / "md_ready.py"),
                            str(cur), "--out", str(dst)], timeout=600)
            if rc == 0 and dst.exists():
                v3, d3 = _verdict(dst)
                if v3 == "PASS":
                    notes.append("清理：修饰残基已降级、水和小分子已去除")
            else:
                shutil.copy2(cur, dst)
                v3 = _verdict(dst)[0]
                notes.append("清理没成功，输出的是上一步的文件")
            _say(f"    成品：{v3} — {dst.name}")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

        rows.append({"file": f.name, "v0": v0, "v1": v1, "v2": v2, "v3": v3,
                     "note": "；".join(notes), "output": dst.name})

    md = ["# 受体准备汇总", "",
          "| 文件 | 原始 | 补 loop | 补侧链 | 清理后 | 说明 | 输出 |",
          "|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f'| {r["file"]} | {r["v0"]} | {r["v1"]} | {r["v2"]} | {r["v3"]} | '
                  f'{r["note"] or "—"} | {r["output"] or "—"} |')
    n = {k: sum(1 for r in rows if r["v3"] == k) for k in ("PASS", "WARN", "FAIL")}
    md += ["", f"共 {len(rows)} 个：PASS {n['PASS']}，WARN {n['WARN']}，FAIL {n['FAIL']}。",
           "", "输出的是 MD 用的干净结构（修饰残基已降级成母体、共晶小分子和水已去掉）。",
           "磷酸化等翻译后修饰在降级时丢失，要保留必须单独参数化。", ""]

    rpt = Path(a.report) if a.report else outdir / "受体准备汇总.md"
    rpt.write_text("\n".join(md), encoding="utf-8")
    _say(f"\n汇总 → {rpt}")
    if a.csv:
        with Path(a.csv).open("w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=["file", "v0", "v1", "v2", "v3", "note", "output"])
            w.writeheader()
            for r in rows:
                w.writerow(r)
        _say(f"CSV → {a.csv}")

    for r in rows:
        _say(f'  {r["v3"]:5s}  {r["file"]}')
    return 0 if n["FAIL"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
