#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 ATB（Automated Topology Builder）网页导出的配体参数整理成流水线能直接用的形状。

ATB 的官方 API 要发邮件申请，网页导出是任何人都能走的路。导出的包长这样：

    1LIG.itp      [ atomtypes ] + [ moleculetype ] + 键角二面角
    1LIG.gro      配体坐标（也可能是 .pdb）
    posre.itp     位置约束（不一定有）

本脚本干四件事，都是 ATB 导出包直接丢给 grompp 会踩的坑：

1. **去掉与蛋白力场重复的 atomtype。** ATB 的 GROMOS 54A7 导出会在 itp 顶上重新
   定义一遍 `C`/`CH2`/`OA` 这些类型，而 gromos54a7.ff 里本来就有，grompp 会报
   atom type 冲突。凡是在目标力场的 ffnonbonded.itp 里已经定义过的类型，这里删掉，
   让 grompp 去力场里取；力场没有的（比如 GAFF2 的 c3/os/ho）原样保留。
2. **统一 moleculetype 名和残基名**，默认改成 LIG，长度压在 5 字符以内（gro 的列宽）。
3. **校对**：itp 的原子数 / 原子顺序 / 原子名要和 gro 一一对上，总电荷要打印出来。
4. **补位置约束**：没给 posre.itp 就按 gro 生成一份 1000 kJ/mol/nm² 的。

用法：
    python scripts/atb_import.py ATB导出目录/ --out 02_topology --forcefield gromos54a7
    python scripts/atb_import.py ligand.zip --check          # 只看诊断，不写文件

输出：`<res>.itp`、`<res>.gro`、`posre_<res>.itp`、`<res>_import.json`
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path

RESMAX = 5          # gro 的残基名列宽
NAMEMAX = 5         # gro 的原子名列宽


# ---------------------------------------------------------------- 读文件
def read_text(p: Path) -> str:
    for enc in ("utf-8-sig", "utf-8", "gbk", "latin-1"):
        try:
            return p.read_text(encoding=enc)
        except UnicodeDecodeError:
            continue
    return p.read_text(encoding="latin-1", errors="replace")


def unpack(src: Path, workdir: Path) -> Path:
    """目录直接用；zip/tgz 解到临时目录。"""
    if src.is_dir():
        return src
    name = src.name.lower()
    if name.endswith(".zip"):
        with zipfile.ZipFile(src) as zf:
            zf.extractall(workdir)
    elif name.endswith((".tgz", ".tar.gz", ".tar.bz2", ".tar")):
        with tarfile.open(src) as tf:
            tf.extractall(workdir, filter="data")
    else:
        raise SystemExit(f"看不懂的输入：{src}（要目录、.zip 或 .tgz）")
    return workdir


# ---------------------------------------------------------------- 解析 itp
def sections(text: str) -> list[tuple[str, list[str]]]:
    """按 [ 段名 ] 切块，返回 [(段名, 数据行)]，注释行去掉。"""
    out: list[tuple[str, list[str]]] = []
    cur, buf = None, []
    for raw in text.splitlines():
        s = raw.split(";", 1)[0].rstrip()
        if not s.strip():
            continue
        m = re.match(r"\s*\[\s*([^\]]+?)\s*\]", s)
        if m:
            if cur is not None:
                out.append((cur, buf))
            cur, buf = m.group(1).strip().lower(), []
            continue
        if cur is not None:
            buf.append(s)
    if cur is not None:
        out.append((cur, buf))
    return out


def atom_types_of(text: str) -> list[str]:
    for name, lines in sections(text):
        if name == "atomtypes":
            return [ln.split()[0] for ln in lines if ln.split()]
    return []


def ff_atomtypes(ffdir: Path) -> set[str]:
    names: set[str] = set()
    for f in ("ffnonbonded.itp", "ffnonbonded_*.itp"):
        for p in sorted(ffdir.glob(f)):
            names |= set(atom_types_of(read_text(p)))
    if not names:                      # 老式力场把 atomtypes 写在 ff*.itp 里
        for p in sorted(ffdir.glob("ff*.itp")):
            names |= set(atom_types_of(read_text(p)))
    return names


def atoms_block(text: str) -> list[list[str]]:
    for name, lines in sections(text):
        if name == "atoms":
            return [ln.split() for ln in lines if ln.split()]
    return []


# ---------------------------------------------------------------- 改写 itp
def rewrite_itp(text: str, res: str, ff_types: set[str]) -> tuple[str, dict]:
    """去掉与力场重复的 atomtypes、改 moleculetype 与残基名。返回 (新文本, 统计)。"""
    secs = sections(text)
    stats = {"dropped": [], "kept": [], "old_molname": None}

    chunks: list[str] = [
        "; 由 atb_import.py 整理：atomtypes 已按目标力场去重，moleculetype 与残基名已统一",
    ]
    for name, lines in secs:
        if name == "atomtypes":
            keep = []
            for ln in lines:
                t = ln.split()[0]
                if t in ff_types:
                    stats["dropped"].append(t)
                else:
                    keep.append(ln)
                    stats["kept"].append(t)
            if not keep:
                continue               # 全部来自力场，整个段删掉
            chunks.append("[ atomtypes ]")
            chunks.extend(keep)
            continue

        if name == "moleculetype":
            data = [ln for ln in lines if ln.split()]
            if data:
                stats["old_molname"] = data[0].split()[0]
                nrexcl = data[0].split()[1] if len(data[0].split()) > 1 else "3"
                chunks.append("[ moleculetype ]")
                chunks.append(f" {res:<8s} {nrexcl}")
                chunks.extend(data[1:])
            else:
                chunks.append("[ moleculetype ]")
                chunks.append(f" {res:<8s} 3")
            continue

        if name == "atoms":
            chunks.append("[ atoms ]")
            for ln in lines:
                q = ln.split()
                if len(q) >= 5:
                    q[3] = res                 # 残基名
                    # 列宽对齐：nr type resnr resname atom ...
                    ln = "".join(f"{q[0]:>6s} {q[1]:<6s} {q[2]:>4s} {q[3]:<5s} "
                                 + " ".join(f"{x:>9s}" for x in q[4:]))
                chunks.append(ln)
            continue

        chunks.append(f"[ {name} ]")
        chunks.extend(lines)
    return "\n".join(chunks) + "\n", stats


def _check_masses(it_atoms: list[list[str]], dropped: list[str]):
    """GROMOS 的原子类型本身质量写 0，真实质量在 [ atoms ] 的第 8 列。
    删掉 itp 的 atomtypes 之后，如果原子行没带质量，原子就会变成无质量——
    这种错误 grompp 不报，跑起来才发现动力学不对。"""
    no_mass_col, zero = [], []
    for q in it_atoms:
        if len(q) < 8:
            no_mass_col.append(q[4])
            continue
        try:
            m = float(q[7])
        except ValueError:
            no_mass_col.append(q[4])
            continue
        if m <= 1e-4 and not q[1].upper().startswith("DUM"):
            zero.append(q[4])
    if no_mass_col:
        print(f"注意：{len(no_mass_col)} 个原子没有质量列（{', '.join(no_mass_col[:5])}），"
              f"质量要靠 atomtype 提供；本次删掉了 {len(set(dropped))} 个力场已有的类型，"
              f"确认它们在 {dropped[0] if dropped else '力场'} 里的质量不是 0。")
    if zero:
        print(f"警告：{len(zero)} 个原子质量为 0（{', '.join(zero[:5])}）——"
              f"grompp 不会报错，但动力学是错的。")


# ---------------------------------------------------------------- 读 gro
def read_gro(p: Path) -> tuple[list[tuple[str, str]], list[str]]:
    """返回 ([(resname, atomname)], 原始行)。"""
    lines = read_text(p).splitlines()
    n = int(lines[1].split()[0])
    body = lines[2:2 + n]
    atoms = [(ln[5:10].strip(), ln[10:15].strip()) for ln in body]
    return atoms, body


def write_gro(src_body: list[str], title: str, res: str, out: Path):
    fixed = []
    for ln in src_body:
        ln = ln.ljust(44)
        fixed.append(f"{ln[:5]}{res:<5s}{ln[10:]}")
    box = "   2.00000   2.00000   2.00000\n"     # 占位盒子，合并时会被蛋白盒子覆盖
    out.write_text(f"{title}\n{len(fixed)}\n" + "\n".join(fixed) + "\n" + box,
                   encoding="utf-8")


def pdb_to_gro(pdb: Path, gro: Path):
    """没有 gro 时用 gmx editconf 转一个（盒子先随便给 2 nm）。"""
    import subprocess
    gmx = shutil.which("gmx") or "gmx"
    r = subprocess.run([gmx, "editconf", "-f", str(pdb), "-o", str(gro),
                        "-c", "-box", "2", "2", "2"],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        raise SystemExit(f"gmx editconf 转换失败：\n{r.stdout[-800:]}\n{r.stderr[-800:]}")


def write_posre(gro: Path, res: str, out: Path):
    atoms, _ = read_gro(gro)
    lines = ["; position restraints generated by atb_import.py", "",
             "[ position_restraints ]", "; atom  type      fx      fy      fz"]
    lines += [f"{i:6d}     1  1000  1000  1000" for i in range(1, len(atoms) + 1)]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------- 找文件
def pick(src: Path) -> tuple[Path, Path | None, Path | None]:
    """在导出包里认出 配体 itp / gro / posre。"""
    itps = [p for p in src.rglob("*.itp")
            if "[ moleculetype ]" in read_text(p).lower().replace("  ", " ")]
    if not itps:
        raise SystemExit(f"{src} 里没找到含 [ moleculetype ] 的 .itp")
    # ATB 的文件名通常是 <名字>.itp；posre 那份不含 moleculetype，已被过滤
    itp = sorted(itps, key=lambda p: (len(p.name), p.name))[0]

    gro = None
    for pat in ("*.gro",):
        cand = sorted(src.rglob(pat))
        if cand:
            gro = cand[0]
    pdb = None
    if gro is None:
        cand = sorted(src.rglob("*.pdb"))
        if cand:
            pdb = cand[0]

    posre = None
    for p in sorted(src.rglob("*.itp")):
        if p == itp:
            continue
        if "[ position_restraints ]" in read_text(p).lower():
            posre = p
            break
    return itp, gro or pdb, posre


# ---------------------------------------------------------------- 主流程
def main(argv=None) -> int:
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    ap = argparse.ArgumentParser(description="整理 ATB 导出的配体参数，供 mdpipe 使用")
    ap.add_argument("src", help="ATB 导出目录，或 .zip / .tgz")
    ap.add_argument("--out", help="输出目录（默认 02_topology）", default="02_topology")
    ap.add_argument("--resname", default="LIG", help="统一后的残基名（默认 LIG，≤5 字符）")
    ap.add_argument("--forcefield", default="gromos54a7",
                    help="目标力场名，用来查 ffnonbonded.itp 里已有哪些 atomtype")
    ap.add_argument("--gmx", default="gmx", help="gmx 可执行文件（只在需要 pdb→gro 时用）")
    ap.add_argument("--check", action="store_true", help="只诊断，不写文件")
    a = ap.parse_args(argv)

    src = Path(a.src).expanduser()
    if not src.exists():
        raise SystemExit(f"找不到 {src}")
    res = a.resname.strip()
    if not (1 <= len(res) <= RESMAX):
        raise SystemExit(f"残基名 {res!r} 长度要在 1–{RESMAX} 之间（gro 列宽限制）")

    tmp = Path(tempfile.mkdtemp(prefix="atb_"))
    try:
        root = unpack(src, tmp)
        itp, coords, posre = pick(root)
        print(f"配体拓扑：{itp.relative_to(root) if root in itp.parents else itp.name}")
        print(f"配体坐标：{coords.name if coords else '（没有）'}")
        print(f"位置约束：{posre.name if posre else '（没有，会生成）'}")

        # 目标力场里的 atomtype 清单
        ff_types: set[str] = set()
        ffdir = None
        gmxp = Path(a.gmx)
        if gmxp.exists():
            for up in (gmxp.parent.parent, gmxp.parent):
                cand = up / "share" / "gromacs" / "top" / f"{a.forcefield}.ff"
                if cand.is_dir():
                    ffdir = cand
                    break
        if ffdir is None:
            print(f"注意：没找到 {a.forcefield}.ff，跳过 atomtype 去重", file=sys.stderr)
        else:
            ff_types = ff_atomtypes(ffdir)
            print(f"力场 {a.forcefield}：{len(ff_types)} 个 atomtype")

        text = read_text(itp)
        new_text, stats = rewrite_itp(text, res, ff_types)
        if stats["dropped"]:
            print(f"从 itp 删掉 {len(stats['dropped'])} 个与力场重复的类型："
                  + ", ".join(sorted(set(stats["dropped"]))))
        if stats["kept"]:
            print(f"保留 {len(stats['kept'])} 个力场没有的类型："
                  + ", ".join(sorted(set(stats["kept"]))))

        # 原子校对
        it_atoms = atoms_block(new_text)
        if not it_atoms:
            raise SystemExit("新 itp 里没有 [ atoms ]，原始文件可能不完整")
        names_itp = [q[4] for q in it_atoms]
        try:
            qsum = sum(float(q[6]) for q in it_atoms)
        except (IndexError, ValueError):
            qsum = float("nan")
        print(f"itp 原子数 {len(it_atoms)}，总电荷 {qsum:+.3f}")
        _check_masses(it_atoms, stats["dropped"])

        # 原子类型必须能在"目标力场 + itp 自带 atomtypes"里找到。ATB 的全原子导出
        # （COpt/HOpt/OAlc/HS14 这些）在原生 gromos54a7.ff 里是没有的，必须配 ATB
        # 发布的 gromos54a7_atb.ff；只看电荷和原子数发现不了这一点，grompp 才会报错。
        if ff_types:
            used = {q[1] for q in it_atoms if len(q) > 1}
            defined = ff_types | set(stats["kept"])
            missing = sorted(used - defined)
            if missing:
                print(f"错误：{len(missing)} 个原子类型在 {a.forcefield}.ff 里没有："
                      f"{', '.join(missing)}")
                print("  ATB 的全原子拓扑要配 ATB 的扩展力场（gromos54a7_atb.ff，"
                      "分子页右侧 Forcefield Files 里下载），原生 54A7 里没有这些类型。")
                return 2
            print(f"原子类型全部在力场里（用到 {len(used)} 个：{', '.join(sorted(used))}）")

        if coords is None:
            raise SystemExit("导出包里没有坐标文件（.gro/.pdb），无法校对原子")

        gro_path = coords
        if coords.suffix.lower() == ".pdb":
            gro_path = tmp / "converted.gro"
            pdb_to_gro(coords, gro_path)
        gro_atoms, _ = read_gro(gro_path)
        names_gro = [n for _, n in gro_atoms]
        print(f"坐标原子数 {len(names_gro)}")
        if len(names_gro) != len(it_atoms):
            raise SystemExit(f"原子数不一致：itp {len(it_atoms)}，坐标 {len(names_gro)}")
        bad = [(i + 1, x, y) for i, (x, y) in enumerate(zip(names_itp, names_gro)) if x != y]
        if bad:
            print(f"注意：{len(bad)} 个原子名对不上，前 3 个：{bad[:3]}")
        else:
            print("原子名与顺序完全一致")

        if a.check:
            print("\n--check：没有写文件。")
            return 0

        out = Path(a.out)
        out.mkdir(parents=True, exist_ok=True)
        (out / f"{res}.itp").write_text(new_text, encoding="utf-8")
        _, body = read_gro(gro_path)
        write_gro(body, f"{res} imported from ATB by atb_import.py", res, out / f"{res}.gro")
        if posre is not None:
            shutil.copy2(posre, out / f"posre_{res}.itp")
        else:
            write_posre(out / f"{res}.gro", res, out / f"posre_{res}.itp")

        info = {"source": str(src), "resname": res, "forcefield": a.forcefield,
                "atoms": len(it_atoms), "charge": round(qsum, 4) if qsum == qsum else None,
                "dropped_atomtypes": sorted(set(stats["dropped"])),
                "kept_atomtypes": sorted(set(stats["kept"])),
                "old_moleculetype": stats["old_molname"]}
        (out / f"{res}_import.json").write_text(
            json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"\n已写出 {out}/{res}.itp、{res}.gro、posre_{res}.itp、{res}_import.json")
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
