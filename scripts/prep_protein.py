#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
prep_protein.py —— 晶体结构预处理: 选 altLoc、去杂原子、保留二硫键记录

为什么需要它:
  pdb2gmx 只认 ATOM 记录和 SSBOND 表头。晶体 PDB 里常见
    ① 同一残基有 A/B 两套坐标 (alternate location), 必须选一套;
    ② 配体/水/离子以 HETATM 形式存在, 跑 apo 模拟时要去掉;
    ③ 二硫键信息藏在 SSBOND 表头里, 若按"只留 ATOM"的方式清洗就会丢。
  最坑的情况: SSBOND 指向的恰好是 B 构象, 统一取 A 会把二硫键拉成 4-5 Å 的
  "断键", pdb2gmx 不报错, 拓扑里却成了两个游离半胱氨酸 —— 模拟结果直接错。

用法:
  python prep_protein.py --in 7G1Q.pdb --out protein.pdb
  python prep_protein.py --in 7G1Q.pdb --out protein.pdb --force-altloc 120:B,127:B
"""
import argparse, collections, sys
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chain", default=None, help="只保留指定链, 如 A")
    ap.add_argument("--keep-het", default="", help="保留的 HETATM 残基名, 逗号分隔, 如 ZN,MG")
    ap.add_argument("--force-altloc", default="",
                    help="强制某残基使用某构象, 如 120:B,127:B (保持二硫键/金属配位)")
    a = ap.parse_args()

    force = {}
    for item in filter(None, a.force_altloc.split(",")):
        res, alt = item.split(":")
        force[res.strip()] = alt.strip().upper()
    keep_het = {x.strip().upper() for x in a.keep_het.split(",") if x.strip()}

    lines = open(a.inp, encoding="utf-8", errors="replace").read().splitlines(True)

    # ---- 1. 统计每个残基各 altLoc 的平均占据率 ----
    occ = collections.defaultdict(list)
    for l in lines:
        if l.startswith("ATOM") and l[16] in "ABCD":
            key = (l[21], l[22:26].strip())
            occ[key + (l[16],)].append(float(l[54:60] or 1.0))
    choice = {}
    for (ch, res, alt), v in occ.items():
        if ch is not None and a.chain and ch != a.chain:
            continue
        choice.setdefault((ch, res), {})[alt] = sum(v) / len(v)

    altloc_of = {}
    for key, d in choice.items():
        if len(d) == 1:
            altloc_of[key] = next(iter(d))
            continue
        if key[1] in force:
            altloc_of[key] = force[key[1]]
        else:
            best = max(d.items(), key=lambda kv: (kv[1], -ord(kv[0])))   # 占据高者优先, 平局取 A
            altloc_of[key] = best[0]

    # ---- 2. 过滤原子 ----
    out, dropped_het, n_alt = [], collections.Counter(), 0
    for l in lines:
        if l.startswith("ATOM"):
            if a.chain and l[21] != a.chain:
                continue
            key = (l[21], l[22:26].strip())
            want = altloc_of.get(key)
            if l[16].strip() and want is not None and l[16] != want:
                n_alt += 1
                continue
            out.append(l[:16] + " " + l[17:])          # 去掉 altLoc 字符; 无 altLoc 的原子(主链)始终保留
        elif l.startswith("HETATM"):
            name = l[17:20].strip().upper()
            if name in keep_het:
                out.append(l)
            else:
                dropped_het[name] += 1

    # ---- 3. 保留 SSBOND 表头, 且只保留两端都存在的 ----
    have = {(l[21], l[22:26].strip()) for l in out if l.startswith("ATOM")}
    for l in lines:
        if not l.startswith("SSBOND"):
            continue
        c1, r1, c2, r2 = l[15], l[17:21].strip(), l[29], l[31:35].strip()
        if (c1, r1) in have and (c2, r2) in have:
            out.insert(0, l)

    out += ["TER\n", "END\n"]
    open(a.out, "w").writelines(out)

    print(f"写出 {a.out}: {sum(1 for l in out if l.startswith(('ATOM','HETATM')))} 个原子")
    print(f"丢弃 altLoc 原子 {n_alt} 个, 丢弃 HETATM: {dict(dropped_het) or '无'}")
    print("altLoc 选择:", {f"{k[1]}": v for k, v in sorted(altloc_of.items(), key=lambda kv: int(kv[0][1]))})

    # ---- 4. 检查二硫键实际距离 ----
    sg = {}
    for l in out:
        if l.startswith("ATOM") and l[12:16].strip() == "SG":
            sg[l[22:26].strip()] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    for l in out:
        if l.startswith("SSBOND"):
            r1, r2 = l[17:21].strip(), l[31:35].strip()
            if r1 in sg and r2 in sg:
                d = float(np.linalg.norm(sg[r1] - sg[r2]))
                flag = "OK" if d < 2.5 else "!! 距离异常, 二硫键没接上"
                print(f"SSBOND CYS{r1}-CYS{r2}: SG-SG = {d:.2f} Å  {flag}")

if __name__ == "__main__":
    main()
