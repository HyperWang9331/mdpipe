#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
merge_complex.py —— 把蛋白 .gro 与配体 .gro 合并成复合物 .gro

为什么不能直接 cat: .gro 第二行记着原子总数, 残基编号必须连续,
最后一行是盒子矢量。直接拼接会让 grompp 报 "number of coordinates
does not match topology"。

用法:
  python merge_complex.py --protein protein.gro --ligand LIG.gro \
      --ligand-resname LIG --out complex.gro
"""
import argparse

def read_gro(path):
    lines = open(path).read().splitlines()
    n = int(lines[1])
    atoms = lines[2:2 + n]
    box = lines[2 + n].split()
    return atoms, box

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein", required=True)
    ap.add_argument("--ligand", required=True)
    ap.add_argument("--ligand-resname", default="LIG")
    ap.add_argument("--out", default="complex.gro")
    a = ap.parse_args()

    patoms, pbox = read_gro(a.protein)
    latoms, _ = read_gro(a.ligand)

    def last_resid(lines):
        return int(lines[-1][0:5])

    out = ["Protein-ligand complex", f"{len(patoms) + len(latoms):5d}"]

    # .gro 列宽: resid(5) resname(5) atomname(5) atomnum(5) x(8) y(8) z(8)
    # 原子顺序决定拓扑对应关系, 原子编号本身不参与解析, 因此原样保留 5:20 段
    for l in patoms:
        out.append(l)

    base = last_resid(patoms)
    for l in latoms:
        resid = base + 1                                   # 配体整体作为一个残基
        out.append(f"{resid % 100000:5d}{a.ligand_resname:<5s}{l[10:20]}{l[20:]}")

    box = pbox if any(float(x) > 0 for x in pbox[:3]) else ["0.00000"] * 3
    out.append("".join(f"{float(x):10.5f}" for x in box[:3]))

    open(a.out, "w").write("\n".join(out) + "\n")
    print(f"写出 {a.out}: {len(patoms)} (蛋白) + {len(latoms)} (配体) = "
          f"{len(patoms) + len(latoms)} 个原子")

if __name__ == "__main__":
    main()
