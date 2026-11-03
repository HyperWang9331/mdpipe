#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
assemble_top.py —— 把配体分子类型挂进 pdb2gmx 生成的蛋白拓扑

做的事:
  1. 在水的 #include 之前插入配体 itp (及其位置约束)
  2. 在 [ molecules ] 里登记配体分子数
每插入一项都会打印, 方便核对顺序。

用法:
  python assemble_top.py --top topol_protein.top --lig-itp LIG.itp \
      --lig-resname LIG --lig-count 1 --out topol.top
"""
import argparse, re

ap = argparse.ArgumentParser()
ap.add_argument("--top", required=True)
ap.add_argument("--lig-itp", default="LIG.itp")
ap.add_argument("--lig-resname", default="LIG")
ap.add_argument("--lig-count", type=int, default=1)
ap.add_argument("--out", default="topol.top")
a = ap.parse_args()

lines = open(a.top).read().splitlines(True)

# ---- 1a. 配体 itp 必须放在第一个 [ moleculetype ] 之前 ----
# GROMACS 要求 [ atomtypes ] 出现在任何 [ moleculetype ] 之前, 否则报
# "Invalid order for directive atomtypes"。配体 itp 里带 [ atomtypes ], 所以
# 紧跟在 forcefield.itp 之后插入。
ff = next(i for i, l in enumerate(lines)
          if l.startswith("#include") and "forcefield.itp" in l)
lines[ff + 1:ff + 1] = ["\n; ---- ligand ----\n", f'#include "{a.lig_itp}"\n']

# ---- 1b. 位置约束放在水的 include 之前 ----
anchor = None
for i, l in enumerate(lines):
    if l.startswith("#include") and re.search(r"(tip3p|tip4p|spc|tip5p)\.itp", l):
        anchor = i
        break
if anchor is None:
    raise SystemExit("找不到水模型的 #include 行, 请检查拓扑文件")

block = [
    f'#ifdef POSRES\n#include "posre_{a.lig_resname}.itp"\n#endif\n',
]
lines[anchor:anchor] = block

# ---- 2. [ molecules ] 登记 ----
mi = next(i for i, l in enumerate(lines) if l.strip() == "[ molecules ]")
end = mi + 1
while end < len(lines) and not lines[end].startswith("["):
    end += 1
# 写在 [ molecules ] 段末尾
insert_at = end
lines.insert(insert_at, f"{a.lig_resname:<16s}{a.lig_count}\n")

open(a.out, "w").writelines(lines)
print(f"写出 {a.out}")
print("  插入:", "".join(block).replace("\n", " ").strip())
print(f"  [ molecules ] 追加: {a.lig_resname} {a.lig_count}")
