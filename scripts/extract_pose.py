#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extract_pose.py —— 从 Vina/AutoDock 的 output.pdbqt 里取出指定构象并补全氢

背景:
  Vina 输出的位姿只有重原子 + 极性氢 (非极性氢被合并掉了)。直接拿去
  MD 会因为缺原子、缺键级而无法参数化; 必须先用 OpenBabel 按化学规则
  补回全部氢原子。

用法:
  python extract_pose.py --pose output.pdbqt --model 1 --out lig_pose.pdb
  python extract_pose.py --pose output.pdbqt --model 1 --out lig_pose.pdb --addh lig_pose.mol2

--addh 会调用 obabel -h 生成带全氢的 mol2, 可直接送 antechamber。
"""
import argparse, subprocess, sys

def read_models(path):
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    blocks, cur = [], None
    for l in lines:
        if l.startswith("MODEL"):
            cur = []
        elif l.startswith("ENDMDL"):
            if cur is not None:
                blocks.append(cur)
            cur = None
        elif cur is not None and l.startswith(("ATOM", "HETATM")):
            cur.append(l)
    if not blocks:
        single = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
        if single:
            blocks = [single]
    return blocks

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pose", required=True, help="Vina 的 output.pdbqt")
    ap.add_argument("--model", type=int, default=1, help="第几个构象 (1 起)")
    ap.add_argument("--resname", default="LIG")
    ap.add_argument("--out", default="lig_pose.pdb")
    ap.add_argument("--addh", default=None, help="补氢后的 mol2 输出文件名")
    a = ap.parse_args()

    blocks = read_models(a.pose)
    if not blocks:
        sys.exit(f"{a.pose} 里没有找到构象")
    if a.model > len(blocks):
        sys.exit(f"文件只有 {len(blocks)} 个构象, 请求的是第 {a.model} 个")

    out = []
    for l in blocks[a.model - 1]:
        out.append(l[:17] + f"{a.resname:<3s}" + l[20:])
    out.append("END")
    open(a.out, "w").write("\n".join(out) + "\n")
    print(f"取出第 {a.model}/{len(blocks)} 个构象 -> {a.out} "
          f"({len(blocks[a.model-1])} 个原子, 仅重原子+极性氢)")

    if a.addh:
        cmd = f'obabel "{a.out}" -O "{a.addh}" -h'
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"obabel 补氢失败:\n{r.stderr}")
        print(f"obabel 补氢 -> {a.addh}")
        print("  下一步请核对: 原子数应为", "重原子+全氢")
        print("  建议再用 RDKit 检查 SMILES 是否与预期分子一致")

if __name__ == "__main__":
    main()
