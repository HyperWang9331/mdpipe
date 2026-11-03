#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""分析结果自检。

只用与 PBC 无关的量（Rg、SASA、接触原子数）去交叉验证 RMSD 类指标，
专门抓那次 `-pbc nojump` 造成的假象：骨架 RMSD 突然大 4 倍，而 Rg/SASA
看着完全正常。自检同时核对时间轴与帧数，防止分析跑到一半或者读错轨迹。

单独用法：
    python scripts/selfcheck.py 07_analysis_r1 [--prod-ns 100]
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np

OK, WARN, FAIL = "ok", "warn", "fail"
RANK = {OK: 0, WARN: 1, FAIL: 2}

# 合理区间（FABP5，~130 残基，单体）
RG_RANGE = (1.30, 1.60)          # nm
SASA_RANGE = (55.0, 115.0)       # nm^2
BACKBONE_TAIL_MAX = 0.40         # nm，尾段均值超过就 warn
BACKBONE_TAIL_FAIL = 0.70        # nm
LIGAND_TAIL_MAX = 1.20           # nm，配体尾段位移
CONTACT_LOW, CONTACT_HIGH = 200, 700
BB_SPIKE_RATIO = 4.0             # 最大值 / 尾段均值


def _xvg(path: pathlib.Path):
    if not path.exists():
        return None
    try:
        d = np.loadtxt(path, comments=["#", "@"])
    except Exception:
        return None
    if d.ndim == 1:
        d = d.reshape(1, -1)
    return d if d.size else None


def _tail(y, frac=0.5):
    k = max(1, int(len(y) * frac))
    return y[-k:]


def check_replicate(ana: pathlib.Path, prod_ns: float | None = None,
                    dt_ns: float = 0.02) -> dict:
    ana = pathlib.Path(ana)
    rep = 1
    m = [p for p in ana.glob("rmsd_backbone_*.xvg")]
    if m:
        try:
            rep = int(m[0].stem.rsplit("_", 1)[1])
        except ValueError:
            rep = 1
    checks: list[dict] = []

    def add(name, verdict, detail):
        checks.append({"check": name, "verdict": verdict, "detail": detail})

    files = {
        "backbone": ana / f"rmsd_backbone_{rep}.xvg",
        "ligand": ana / f"rmsd_ligand_{rep}.xvg",
        "rg": ana / f"rg_{rep}.xvg",
        "sasa": ana / f"sasa_{rep}.xvg",
        "contacts": ana / "contacts.xvg",
        "hbond": ana / f"hbond_{rep}.xvg",
    }
    missing = [k for k, p in files.items() if not p.exists()]
    if missing:
        add("文件齐全", FAIL, "缺 " + "、".join(missing))
        return _finish(checks, rep)

    data = {k: _xvg(p) for k, p in files.items()}
    t = data["backbone"][:, 0]

    # 1. 时间轴
    dt = np.diff(t)
    if np.any(dt <= 0):
        add("时间轴单调", FAIL, "存在非递增的时间点")
    else:
        add("时间轴单调", OK, f"{t[0]:.2f} → {t[-1]:.2f} ns")
    # xvg 里的时间是 %g 格式，0.02 这种值会有末位抖动，容差按 1 ps 或 2% 取大
    tol = max(1e-3, 0.02 * float(dt.mean())) if dt.size else 1e-3
    if dt.size and (dt.max() - dt.min()) > tol:
        add("时间间隔均匀", WARN, f"dt 在 {dt.min():.3f}–{dt.max():.3f} ns 之间变化")
    else:
        add("时间间隔均匀", OK, f"每帧 {dt.mean():.3f} ns" if dt.size else "单帧")
    if prod_ns and abs(t[-1] - prod_ns) > max(0.02 * prod_ns, 0.1):
        add("轨迹长度", WARN, f"末帧 {t[-1]:.2f} ns，任务是 {prod_ns:g} ns")
    else:
        add("轨迹长度", OK, f"{t[-1]:.2f} ns / {len(t)} 帧")

    # 2. 首帧与尾段骨架 RMSD
    bb = data["backbone"][:, 1]
    if bb[0] > 0.05:
        add("首帧基准", WARN, f"首帧骨架 RMSD = {bb[0]:.3f} nm（应≈0）")
    else:
        add("首帧基准", OK, f"{bb[0]:.3f} nm")
    bb_tail = _tail(bb)
    bb_mean = float(bb_tail.mean())
    if bb_mean > BACKBONE_TAIL_FAIL:
        add("骨架稳定", FAIL, f"尾段 {bb_mean:.3f} nm > {BACKBONE_TAIL_FAIL} nm")
    elif bb_mean > BACKBONE_TAIL_MAX:
        add("骨架稳定", WARN, f"尾段 {bb_mean:.3f} nm 偏大")
    else:
        add("骨架稳定", OK, f"尾段 {bb_mean:.3f} ± {bb_tail.std():.3f} nm")
    if bb.max() > max(BB_SPIKE_RATIO * bb_mean, 0.5):
        add("RMSD 无尖峰", WARN,
            f"最大 {bb.max():.2f} nm 是尾段均值的 {bb.max() / max(bb_mean, 1e-9):.1f} 倍，"
            f"先怀疑镜像/跳变")

    # 3. 与 PBC 无关的量：Rg / SASA
    for key, rng, label, unit in (("rg", RG_RANGE, "回转半径", "nm"),
                                  ("sasa", SASA_RANGE, "SASA", "nm²")):
        y = data[key][:, 1]
        v = float(_tail(y).mean())
        if not (rng[0] <= v <= rng[1]):
            add(label, WARN, f"尾段 {v:.2f} {unit}，超出常见区间 {rng[0]}–{rng[1]}")
        else:
            add(label, OK, f"尾段 {v:.2f} {unit}")

    # 4. 配体 RMSD 与接触原子数必须互相印证
    lig = float(_tail(data["ligand"][:, 1]).mean())
    con = float(_tail(data["contacts"][:, 1]).mean())
    if lig < 1.0 and con < CONTACT_LOW:
        add("配体指标自洽", FAIL,
            f"配体 RMSD {lig:.2f} nm 说它还在原位，接触原子只有 {con:.0f} 个说不在了")
    elif lig > 2.0 and con > CONTACT_HIGH:
        add("配体指标自洽", WARN,
            f"配体 RMSD {lig:.2f} nm 很大但接触原子仍有 {con:.0f} 个，检查参考结构")
    else:
        add("配体指标自洽", OK, f"配体 RMSD {lig:.2f} nm / 接触原子 {con:.0f}")
    if lig > LIGAND_TAIL_MAX:
        add("配体位移", WARN, f"尾段 {lig:.2f} nm，已离开初始位姿")

    return _finish(checks, rep, {"ligand_tail_nm": round(lig, 3),
                                 "contacts_tail": round(con, 1),
                                 "backbone_tail_nm": round(bb_mean, 3)})


def _finish(checks, rep, metrics=None):
    verdict = OK
    for c in checks:
        if RANK[c["verdict"]] > RANK[verdict]:
            verdict = c["verdict"]
    return {"replicate": rep, "verdict": verdict, "checks": checks,
            "metrics": metrics or {}, "passed": verdict == OK}


def worst(*verdicts):
    out = OK
    for v in verdicts:
        if RANK.get(v, 0) > RANK[out]:
            out = v
    return out


def job_verdict(job_dir) -> str | None:
    """读一个作业所有重复的 selfcheck.json，返回最差的一档。"""
    res = []
    for p in sorted(pathlib.Path(job_dir).glob("07_analysis*/selfcheck.json")):
        try:
            res.append(json.loads(p.read_text(encoding="utf-8")).get("verdict", OK))
        except Exception:
            pass
    return worst(*res) if res else None


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description="分析结果自检")
    ap.add_argument("dir", help="分析目录，如 07_analysis_r1")
    ap.add_argument("--prod-ns", type=float, default=None)
    ap.add_argument("--write", action="store_true", help="写回 selfcheck.json")
    a = ap.parse_args(argv)
    d = pathlib.Path(a.dir)
    r = check_replicate(d, a.prod_ns)
    print(f"复检结果：{r['verdict']}")
    for c in r["checks"]:
        print(f"  [{c['verdict']:4s}] {c['check']}：{c['detail']}")
    if a.write:
        (d / "selfcheck.json").write_text(json.dumps(r, ensure_ascii=False, indent=2),
                                          encoding="utf-8")
    return 0 if r["verdict"] != FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
