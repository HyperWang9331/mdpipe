#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""核心纯函数的单元测试（不需要 GROMACS、不需要跑 MD）。

    python -m unittest discover -s tests -v

覆盖的是那些"错了也能跑完、只是结果没意义"的地方：周期边界下的接触/
氢键判定、自检的阈值、MM-PBSA 的取样窗口、归档前的上游检查、凭据是否
真的没落成明文。
"""

from __future__ import annotations

import json
import os
import pathlib
import shutil
import sys
import tempfile
import unittest

import numpy as np

SCRIPTS = pathlib.Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))

import autoconfig               # noqa: E402
import doctor                   # noqa: E402
import install_skill             # noqa: E402
import interactions as it       # noqa: E402
import mmpbsa                   # noqa: E402
import mdpipe                   # noqa: E402
import mdpipe_gui               # noqa: E402
import mdview                   # noqa: E402
import notify                   # noqa: E402
import selfcheck                # noqa: E402


class TmpDir(unittest.TestCase):
    def setUp(self):
        self.tmp = pathlib.Path(tempfile.mkdtemp(prefix="mdpipe_test_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write(self, rel: str, text: str) -> pathlib.Path:
        p = self.tmp / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
        return p


# ---------------------------------------------------------------- 结构读取
class TestGroAndMasks(TmpDir):
    def test_read_gro_fixed_columns(self):
        def line(resnr, resname, name):
            return (f"{resnr:5d}{resname:<5s}{name:>5s}{1:5d}"
                    "   0.100   0.200   0.300")

        gro = self.write("t.gro", "\n".join([
            "test system",
            "    3",
            line(1, "LIG", "C1"),
            line(1, "LIG", "O1"),
            line(10, "THR", "OG1"),
            "   1.0 1.0 1.0",
        ]) + "\n")
        atoms = it.read_gro(gro)
        self.assertEqual(atoms, [(1, "LIG", "C1"), (1, "LIG", "O1"),
                                 (10, "THR", "OG1")])

    def test_element_names(self):
        for name, elem in [("OG1", "O"), ("NZ", "N"), ("CA", "C"),
                           ("HG1", "H"), ("C12", "C"), ("", "?")]:
            self.assertEqual(it._element(name), elem, name)

    def test_topology_masks(self):
        atoms = [(1, "LIG", "C1"), (1, "LIG", "O1"), (1, "LIG", "H1"),
                 (2, "ALA", "CA"), (2, "ALA", "N")]
        m = it.topology_masks(atoms)
        self.assertEqual(list(m["lig_heavy"]), [0, 1])
        self.assertEqual(list(m["prot_heavy"]), [3, 4])
        self.assertEqual(list(m["lig_o"]), [1])
        self.assertEqual(list(m["prot_o"]), [4])
        self.assertEqual(list(m["lig_all"]), [0, 1, 2])

    def test_donor_parent_assignment(self):
        # 只有一个 O–H 和一个 C–H；只有前者算供体
        atoms = [(1, "LIG", "O1"), (1, "LIG", "H1"),
                 (1, "LIG", "C1"), (1, "LIG", "H2")]
        xyz = np.array([[[0.0, 0, 0], [0.10, 0, 0], [0.5, 0, 0], [0.60, 0, 0]]])
        masks = it.topology_masks(atoms)
        pairs = it._donors(xyz, atoms, masks, "lig_all")
        self.assertEqual(pairs, [(1, 0)])


# ---------------------------------------------------------------- 接触
class TestContacts(unittest.TestCase):
    def test_occupancy_and_bound_fraction(self):
        atoms = [(10, "ALA", "N"), (10, "ALA", "CA"),
                 (1, "LIG", "C1"), (1, "LIG", "O1")]
        masks = it.topology_masks(atoms)
        f1 = np.array([[0, 0, 0], [0, 0, 0.10], [0, 0, 0.40], [0.1, 0, 0.40]])
        f2 = f1.copy()
        f2[2] += [0, 0, 5.0]                      # 配体走远
        f2[3] += [0, 0, 5.0]
        rows, bound = it.contact_table(np.array([f1, f2]), atoms, masks,
                                       cutoff=0.45)
        self.assertEqual(list(bound), [True, False])
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertAlmostEqual(r["occ"], 0.5)
        self.assertAlmostEqual(r["occ_bound"], 1.0)   # 只在结合帧里算
        self.assertAlmostEqual(r["min_dist"], 0.30, places=6)
        self.assertEqual(r["lig_atom"], "C1")

    def test_no_contact_gives_zero_bound(self):
        atoms = [(10, "ALA", "N"), (1, "LIG", "C1")]
        masks = it.topology_masks(atoms)
        xyz = np.array([[[0, 0, 0], [1.0, 0, 0]]])
        rows, bound = it.contact_table(xyz, atoms, masks)
        self.assertFalse(bound.any())
        self.assertEqual(rows[0]["occ"], 0.0)
        self.assertEqual(rows[0]["occ_bound"], 0.0)   # 不除以 0


# ---------------------------------------------------------------- 氢键
def hb_system(acceptor_xyz=(0.30, 0.0, 0.0)):
    atoms = [(10, "THR", "OG1"), (10, "THR", "HG1"),
             (1, "LIG", "O1")]
    masks = it.topology_masks(atoms)
    xyz = np.array([[[0.0, 0.0, 0.0], [0.10, 0.0, 0.0], list(acceptor_xyz)]])
    return atoms, masks, xyz


class TestHbonds(unittest.TestCase):
    def test_linear_hbond_counted(self):
        atoms, masks, xyz = hb_system()
        rows = it.hbond_table(xyz, atoms, masks)
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual((r["resnr"], r["resname"], r["prot_atom"], r["lig_atom"]),
                         (10, "THR", "OG1", "O1"))
        self.assertAlmostEqual(r["occ"], 1.0)
        self.assertAlmostEqual(r["dist"], 0.30, places=6)
        self.assertGreater(r["angle"], 179.0)

    def test_angle_criterion_rejects(self):
        # D···A = 0.316 nm 满足距离，但 D–H···A ≈ 90°，不是氢键
        atoms, masks, xyz = hb_system((0.10, -0.30, 0.0))
        self.assertEqual(it.hbond_table(xyz, atoms, masks), [])

    def test_distance_criterion_rejects(self):
        atoms, masks, xyz = hb_system((0.40, 0.0, 0.0))
        self.assertEqual(it.hbond_table(xyz, atoms, masks), [])

    def test_bound_occupancy_uses_bound_frames(self):
        atoms, masks, xyz = hb_system()
        far = xyz.copy()
        far[0, 2] = [5.0, 5.0, 5.0]
        rows = it.hbond_table(np.array([xyz[0], far[0]]), atoms, masks,
                              bound=np.array([True, False]))
        self.assertAlmostEqual(rows[0]["occ"], 0.5)
        self.assertAlmostEqual(rows[0]["occ_bound"], 1.0)


# ---------------------------------------------------------------- 盐桥
class TestSaltBridges(unittest.TestCase):
    def system(self):
        atoms = [(1, "LIG", "C1"), (1, "LIG", "O1"), (1, "LIG", "O2"),
                 (1, "LIG", "O3"),
                 (50, "ARG", "NE"), (50, "ARG", "NH2")]
        masks = it.topology_masks(atoms)
        xyz = np.array([[[0, 0, 0], [0.0, 0.25, 0], [0.0, -0.25, 0],
                         [1.2, 0, 0],
                         [0.0, 0.25, 0.30], [0.0, -0.25, 0.30]]])
        return atoms, masks, xyz

    def test_carboxylate_only(self):
        atoms, masks, xyz = self.system()
        rows = it.salt_bridges(xyz, atoms, masks, neg_names={"O1", "O2"})
        self.assertTrue(rows)
        self.assertEqual({r["lig_atom"] for r in rows}, {"O1", "O2"})
        self.assertTrue(all(r["resname"] == "ARG" for r in rows))

    def test_without_anionic_set_returns_nothing(self):
        atoms, masks, xyz = self.system()
        self.assertEqual(it.salt_bridges(xyz, atoms, masks, neg_names=set()), [])

    def test_far_apart_not_counted(self):
        atoms, masks, xyz = self.system()
        xyz[0, 4] += [0, 0, 1.0]
        xyz[0, 5] += [0, 0, 1.0]
        self.assertEqual(it.salt_bridges(xyz, atoms, masks,
                                         neg_names={"O1", "O2"}), [])

    def test_ligand_anionic_atoms(self):
        try:
            from rdkit import Chem
        except Exception:
            self.skipTest("没有 RDKit")
        glucuronic = Chem.MolFromSmiles("O=C([O-])C1OC(O)C(O)C(O)C1O")
        glu = Chem.MolFromSmiles("OC1OC(O)C(O)C(O)C1O")
        got = it.ligand_anionic_atoms((glucuronic, None, None, None))
        self.assertEqual(len(got), 2)
        self.assertEqual(it.ligand_anionic_atoms((glu, None, None, None)), set())

    def test_map_ligand_names_rejects_mismatch(self):
        try:
            from rdkit import Chem
        except Exception:
            self.skipTest("没有 RDKit")
        mol = Chem.RemoveHs(Chem.MolFromSmiles("OCC"))
        lig = (mol, np.zeros((3, 2)), {0: 0, 1: 1, 2: 2}, ["O", "C", "C"])
        atoms = [(1, "LIG", "O1"), (1, "LIG", "C1"), (1, "LIG", "C2")]
        self.assertEqual(it.map_ligand_names(atoms, "LIG", lig),
                         {"O1": 0, "C1": 1, "C2": 2})
        bad = [(1, "LIG", "O1"), (1, "LIG", "C1")]      # 少了原子
        self.assertEqual(it.map_ligand_names(bad, "LIG", lig), {})


# ---------------------------------------------------------------- MM-PBSA
class TestMMPBSA(TmpDir):
    def test_win_to_wsl(self):
        self.assertEqual(mmpbsa.win_to_wsl(r"C:\Users\x\job"),
                         "/mnt/c/Users/x/job")
        self.assertEqual(mmpbsa.win_to_wsl("/mnt/e/MD"), "/mnt/e/MD")

    def test_ndx_groups_and_index(self):
        ndx = self.write("index.ndx", "\n".join([
            "[ System ]", "1 2 3", "[ Protein ]", "1 2", "[ LIG ]", "3", ""]))
        self.assertEqual(mmpbsa.ndx_groups(ndx), ["System", "Protein", "LIG"])
        self.assertEqual(mmpbsa.group_index(ndx, "LIG"), 2)
        with self.assertRaises(RuntimeError):
            mmpbsa.group_index(ndx, "Water")

    def test_default_window_is_last_20ns(self):
        # 5001 帧、100 ns：20 ns 窗口应从第 4001 帧开始，每 100 帧取一帧
        out = self.tmp
        mmpbsa.write_input(out, frames=10, nframe=5001, last_ps=100000.0,
                           last_ns=20.0, istrng=0.15, decomp=False)
        text = (out / "mmpbsa.in").read_text(encoding="utf-8")
        self.assertIn("startframe=4001", text)
        self.assertIn("endframe=5001", text)
        self.assertIn("interval=100", text)
        self.assertNotIn("&decomp", text)

    def test_full_trajectory_window(self):
        mmpbsa.write_input(self.tmp, frames=10, nframe=5001, last_ps=100000.0,
                           last_ns=0, istrng=0.15, decomp=False)
        text = (self.tmp / "mmpbsa.in").read_text(encoding="utf-8")
        self.assertIn("startframe=1", text)
        self.assertIn("interval=500", text)

    def test_short_last_window_never_empty(self):
        # 要求 200 ns 但轨迹只有 100 ns：窗口退化为全程，不能出现 start>end
        mmpbsa.write_input(self.tmp, frames=10, nframe=5001, last_ps=100000.0,
                           last_ns=200.0, istrng=0.15, decomp=False)
        text = (self.tmp / "mmpbsa.in").read_text(encoding="utf-8")
        self.assertIn("startframe=1", text)

    def test_env_paths_from_config(self):
        orig = mmpbsa._pipe_config
        try:
            mmpbsa._pipe_config = lambda: {
                "gmx": "E:/gmx/bin/gmx.exe", "mmpbsa_env": "/opt/mmpbsa",
                "wsl_distro": "Debian", "wsl_home": "/home/u"}
            e = mmpbsa.env_paths()
            self.assertEqual(e["env"], "/opt/mmpbsa")
            self.assertEqual(e["distro"], "Debian")
            self.assertEqual(e["home"], "/home/u")
            # 没给 gmx_share 时从 gmx 的安装位置推
            self.assertEqual(e["gmxlib"], "/mnt/e/gmx/share/gromacs/top")

            mmpbsa._pipe_config = lambda: {"gmx_share": "/opt/top"}
            self.assertEqual(mmpbsa.env_paths()["gmxlib"], "/opt/top")
        finally:
            mmpbsa._pipe_config = orig

    def test_parse_results(self):
        dat = self.write("FINAL_RESULTS.dat", "\n".join([
            "|something|",
            "POISSON BOLTZMANN:",
            "-" * 40,
            # 与 gmx_MMPBSA 实测输出一致：空格分隔、无等号、带 SD 列
            "ΔVDWAALS                -15.60          0.10       2.08",
            "ΔEEL                    -49.86          0.30       3.10",
            "ΔEPB                     70.67          0.40       4.00",
            "ΔENPOLAR                 -2.56          0.01       0.10",
            "ΔGGAS                   -65.46          0.50       5.00",
            "ΔGSOLV                   68.11          0.39       3.90",
            "ΔTOTAL                    2.64          2.68       6.01",
            "",
        ]) + "\n")
        r = mmpbsa.parse_results(dat)
        self.assertAlmostEqual(r["dG"], 2.64)
        self.assertAlmostEqual(r["sd"], 2.68)
        self.assertAlmostEqual(r["eel"], -49.86)
        self.assertAlmostEqual(r["epb"], 70.67)
        self.assertEqual(mmpbsa.parse_results(self.tmp / "nope.dat"), {})


# ---------------------------------------------------------------- 自检
def xvg(path: pathlib.Path, t, y, title="x", ylab="y"):
    lines = [f"# {title}", f'@    title "{title}"',
             f'@    yaxis  label "{ylab}"']
    lines += [f"{a:.4f} {b:.4f}" for a, b in zip(t, y)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class TestSelfcheck(TmpDir):
    def make_analysis(self, bb):
        d = self.tmp / "07_analysis_r1"
        d.mkdir(parents=True)
        t = np.linspace(0, 100, 1001)
        xvg(d / "rmsd_backbone_1.xvg", t, bb(t))
        xvg(d / "rmsd_ligand_1.xvg", t, np.full_like(t, 0.2))
        xvg(d / "rg_1.xvg", t, np.full_like(t, 1.45))
        xvg(d / "sasa_1.xvg", t, np.full_like(t, 80.0))
        xvg(d / "contacts.xvg", t, np.full_like(t, 400.0))
        xvg(d / "hbond_1.xvg", t, np.full_like(t, 1.5))
        return d

    def test_stable_trajectory_passes(self):
        d = self.make_analysis(lambda t: 0.02 + 0.0005 * t)
        r = selfcheck.check_replicate(d, prod_ns=100)
        self.assertEqual(r["verdict"], selfcheck.OK, r)
        self.assertEqual(r["replicate"], 1)

    def test_backbone_blowup_fails(self):
        d = self.make_analysis(lambda t: 0.02 + 0.02 * t)   # 末段 2 nm
        r = selfcheck.check_replicate(d, prod_ns=100)
        self.assertEqual(r["verdict"], selfcheck.FAIL)
        names = {c["check"] for c in r["checks"] if c["verdict"] == selfcheck.FAIL}
        self.assertIn("骨架稳定", names)

    def test_missing_files_fail(self):
        d = self.make_analysis(lambda t: np.full_like(t, 0.02))
        (d / "rg_1.xvg").unlink()
        r = selfcheck.check_replicate(d)
        self.assertEqual(r["verdict"], selfcheck.FAIL)
        self.assertEqual(r["checks"][0]["check"], "文件齐全")

    def test_ligand_metric_consistency(self):
        d = self.make_analysis(lambda t: np.full_like(t, 0.02))
        t = np.linspace(0, 100, 1001)
        xvg(d / "rmsd_ligand_1.xvg", t, np.full_like(t, 0.2))
        xvg(d / "contacts.xvg", t, np.full_like(t, 20.0))
        r = selfcheck.check_replicate(d)
        self.assertEqual(r["verdict"], selfcheck.FAIL)

    def test_worst_and_job_verdict(self):
        self.assertEqual(selfcheck.worst("ok", "warn", "fail"), "fail")
        self.assertEqual(selfcheck.worst("ok", "warn"), "warn")
        for rep, verdict in ((1, "ok"), (2, "warn")):
            p = self.tmp / f"07_analysis_r{rep}" / "selfcheck.json"
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(json.dumps({"verdict": verdict}), encoding="utf-8")
        self.assertEqual(selfcheck.job_verdict(self.tmp), "warn")
        self.assertIsNone(selfcheck.job_verdict(self.tmp / "empty"))


# ---------------------------------------------------------------- 归档
class TestArchive(TmpDir):
    def make_job(self):
        j = self.tmp / "jobs" / "J1"
        for rel in ("06_prod/md.xtc", "06_prod/md.tpr",
                    "07_analysis/md_center.xtc", "07_analysis/md_prot.xtc"):
            p = j / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(b"x")
        return j

    def test_complete_production_allows_deletion(self):
        j = self.make_job()
        self.assertEqual(mdpipe._archive_blockers(j), {})
        plan = mdpipe._archive_plan(["J1"], base=self.tmp / "jobs")
        self.assertEqual(len(plan), 2)
        self.assertTrue(all(p["blocked"] is None for p in plan))

    def test_missing_tpr_blocks(self):
        j = self.make_job()
        (j / "06_prod" / "md.tpr").unlink()
        blocked = mdpipe._archive_blockers(j)
        self.assertIn(j / "07_analysis" / "md_center.xtc", blocked)
        self.assertIn("md.tpr", blocked[j / "07_analysis" / "md_center.xtc"])

    def test_unfinished_equilibration_blocks(self):
        j = self.make_job()
        p = j / "04_nvt" / "nvt.xtc"
        p.parent.mkdir(parents=True)
        p.write_bytes(b"x")
        blocked = mdpipe._archive_blockers(j)
        self.assertIn("npt.gro", blocked[p])


# ---------------------------------------------------------------- 安装清单
class TestInstallPlan(unittest.TestCase):
    def test_requirements_table_valid(self):
        levels = {"必需", "默认路线", "建议", "可选"}
        modes = {"自动", "半自动", "手动"}
        for name, level, mode, how, _note in doctor.REQUIREMENTS:
            self.assertIn(level, levels, name)
            self.assertIn(mode, modes, name)
            self.assertTrue(str(how).strip(), name)

    def test_plan_text_has_all_sections(self):
        text = doctor._plan_text()
        for kw in ("可以自动装", "半自动", "只能手动", "最小可跑组合"):
            self.assertIn(kw, text)

    def test_pip_command_shape(self):
        cmd = doctor.pip_command(["numpy", "scipy"])
        self.assertEqual(cmd[1:4], ["-m", "pip", "install"])
        self.assertEqual(cmd[4:], ["numpy", "scipy"])

    def test_missing_packages(self):
        self.assertNotIn("numpy", doctor.missing_packages(["numpy"]))
        bogus = "definitely-not-a-real-package-xyz"
        self.assertIn(bogus, doctor.missing_packages([bogus]))

    def test_cpu_fallback_mdrun(self):
        cfg = {"gmx": "gmx", "ntomp": "4", "gpu_id": 0}
        gpu = mdpipe._mdrun_prefix({**cfg, "use_gpu": True})
        self.assertIn("-nb", gpu)
        self.assertIn("-gpu_id", gpu)
        cpu = mdpipe._mdrun_prefix({**cfg, "use_gpu": False})
        self.assertNotIn("-nb", cpu)
        self.assertNotIn("-gpu_id", cpu)
        self.assertEqual(cpu, ["gmx", "mdrun", "-ntmpi", "1", "-ntomp", "4"])


# ---------------------------------------------------------------- 启动脚本
class TestLaunchers(unittest.TestCase):
    """cmd 解析 .bat 时要求 CRLF；LF 会把中文 rem 行和下一行粘成一条命令。"""

    def test_find_doc_checks_docs_dir(self):
        """安装版文档在 docs\\ 下，菜单不能只找根目录。"""
        with tempfile.TemporaryDirectory() as d:
            root = pathlib.Path(d)
            (root / "docs").mkdir()
            (root / "docs" / "手册.md").write_text("x", encoding="utf-8")
            self.assertEqual(mdpipe_gui._find_doc("手册.md", root),
                             root / "docs" / "手册.md")

    def test_find_doc_missing_returns_none(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertIsNone(mdpipe_gui._find_doc("没有.md", pathlib.Path(d)))

    def test_bat_files_use_crlf(self):
        root = pathlib.Path(__file__).resolve().parent.parent
        bats = sorted(root.glob("installer/*.bat"))
        self.assertTrue(bats, "找不到 installer/*.bat")
        for p in bats:
            raw = p.read_bytes()
            self.assertNotIn(b"\n", raw.replace(b"\r\n", b""),
                             f"{p.name} 里还有裸 LF 换行")
            self.assertIn(b"\r\n", raw, f"{p.name} 没有 CRLF 换行")

    def test_bat_has_no_bom(self):
        root = pathlib.Path(__file__).resolve().parent.parent
        for p in sorted(root.glob("installer/*.bat")):
            self.assertFalse(p.read_bytes().startswith(b"\xef\xbb\xbf"),
                             f"{p.name} 带 BOM 时 cmd 会把首行当命令执行")


# ---------------------------------------------------------------- AI 技能文件
class TestSkill(unittest.TestCase):
    def setUp(self):
        self.root = pathlib.Path(__file__).resolve().parent.parent
        self.skill = self.root / "skills" / "md-pipeline" / "SKILL.md"

    def test_frontmatter(self):
        text = self.skill.read_text(encoding="utf-8")
        self.assertTrue(text.startswith("---\n"))
        front = text.split("---", 2)[1]
        self.assertIn("name: md-pipeline", front)
        self.assertIn("description:", front)

    def test_covers_core_commands(self):
        text = self.skill.read_text(encoding="utf-8")
        for cmd in ("mdpipe doctor", "mdpipe add", "mdpipe run",
                    "mdpipe status", "mdpipe analysis", "mdpipe mmpbsa",
                    "mdpipe archive"):
            self.assertIn(cmd, text)

    def test_install_one_lifecycle(self):
        with tempfile.TemporaryDirectory() as d:
            dest = pathlib.Path(d) / "skills" / "md-pipeline"
            self.assertEqual(install_skill.install_one(dest, "v1"), "installed")
            self.assertEqual(install_skill.install_one(dest, "v1"), "unchanged")
            self.assertEqual(install_skill.install_one(dest, "v2"), "exists")
            self.assertEqual(install_skill.install_one(dest, "v2", force=True),
                             "installed")
            self.assertEqual((dest / "md-pipeline" / "SKILL.md")
                             .read_text(encoding="utf-8"), "v2")

    def test_install_agents_is_idempotent(self):
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d) / "AGENTS.md"
            self.assertEqual(install_skill.install_agents(p), "installed")
            self.assertEqual(install_skill.install_agents(p), "unchanged")
            text = p.read_text(encoding="utf-8")
            self.assertEqual(text.count(install_skill.AGENTS_BEGIN), 1)

    def test_install_agents_keeps_existing_content(self):
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d) / "AGENTS.md"
            p.write_text("原有内容\n", encoding="utf-8")
            self.assertEqual(install_skill.install_agents(p), "exists")
            self.assertEqual(install_skill.install_agents(p, force=True), "installed")
            self.assertTrue(p.read_text(encoding="utf-8").startswith("原有内容"))

    def test_targets_layout(self):
        with tempfile.TemporaryDirectory() as d:
            home = pathlib.Path(d)
            t = install_skill.registry(home=home, cwd=home)
            for key, rel in [
                    ("clawsgo", ".clawsgo/skills"),
                    ("dsh", ".dsh/skills"),
                    ("workbuddy", ".workbuddy/skills"),
                    ("codex", ".codex/skills"),
                    ("claude", ".claude/skills")]:
                self.assertEqual(t[key].path, home.joinpath(*rel.split("/")))
                self.assertEqual(t[key].kind, "skill")
            self.assertEqual(t["trae"].kind, "rule")
            self.assertEqual(t["trae"].path, home / ".trae" / "rules" / "md-pipeline.md")

    def test_install_rule_cursor_frontmatter(self):
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d) / ".cursor" / "rules" / "md-pipeline.mdc"
            body = install_skill.rule_body()
            self.assertEqual(
                install_skill.install_rule(p, body, cursor=True), "installed")
            text = p.read_text(encoding="utf-8")
            self.assertTrue(text.startswith("---\ndescription:"))
            self.assertIn(install_skill.AGENTS_BEGIN, text)
            # 重复写不会叠加 frontmatter
            self.assertEqual(
                install_skill.install_rule(p, body, cursor=True), "unchanged")
            self.assertEqual(text.count("description:"), 1)

    def test_installer_ships_skill(self):
        iss = (self.root / "installer" / "mdpipe.iss").read_text(encoding="utf-8")
        self.assertIn(r"skills\md-pipeline\SKILL.md", iss)
        self.assertIn(r".clawsgo\skills\md-pipeline", iss)


# ---------------------------------------------------------------- Markdown 渲染
class TestMarkdown(unittest.TestCase):
    def test_inline_tokens(self):
        toks = mdview.inline_tokens("普通**粗体**和*斜体*与`code`")
        kinds = [k for k, _, _ in toks]
        self.assertEqual(kinds, ["text", "bold", "text", "italic", "text", "code"])

    def test_link_and_image_tokens(self):
        toks = mdview.inline_tokens("[手册](a.md) ![图](b.png)")
        self.assertEqual(toks[0], ("link", "手册", "a.md"))
        self.assertEqual(toks[2], ("image", "图", "b.png"))

    def test_plain_strips_marks(self):
        self.assertEqual(mdview.plain("**能自动装** 与 [链接](x.md)"), "能自动装 与 链接")

    def test_display_width_cjk(self):
        self.assertEqual(mdview.display_width("ab中"), 4)

    def test_wrap_cell_prefers_punctuation(self):
        lines = mdview.wrap_cell("Ubuntu、Miniforge、AmberTools、Ghostscript", 20)
        self.assertTrue(all(mdview.display_width(x) <= 20 for x in lines))
        self.assertNotIn("G", lines)

    def test_fit_widths_shrinks_to_total(self):
        w = mdview.fit_widths([14, 104, 72], total=92, floor=8)
        self.assertLessEqual(sum(w) + 2 * (len(w) - 1), 92)
        self.assertTrue(all(x >= 8 for x in w))

    def test_split_row(self):
        self.assertEqual(mdview.split_row("| a | b\\|c |"), ["a", "b|c"])


# ---------------------------------------------------------------- 文档检索
class TestDocLookup(unittest.TestCase):
    def test_all_docs_found(self):
        for key in mdpipe.DOC_FILES:
            p = mdpipe.doc_path(key)
            self.assertTrue(p.is_file(), f"{key} -> {p}")

    def test_section_extract(self):
        lines = ["# 标题", "## 1. 甲", "正文甲", "### 1.1 细", "细正文",
                 "## 2. 乙", "正文乙"]
        got = mdpipe._doc_section(lines, "1")
        self.assertEqual(got[0], "## 1. 甲")
        self.assertIn("### 1.1 细", got)
        self.assertNotIn("## 2. 乙", got)

    def test_subsection_extract(self):
        lines = ["## 20. 批", "### 20.10 相互作用", "内容", "### 20.11 换机器", "x"]
        got = mdpipe._doc_section(lines, "20.10")
        self.assertEqual(got, ["### 20.10 相互作用", "内容"])


# ---------------------------------------------------------------- 路径自动探测
class TestAutoConfig(unittest.TestCase):
    def test_pick_distro_prefers_ubuntu(self):
        self.assertEqual(
            autoconfig.pick_distro(["Debian", "Ubuntu-24.04", "docker-desktop"]),
            "Ubuntu-24.04")

    def test_pick_distro_falls_back_to_first(self):
        self.assertEqual(autoconfig.pick_distro(["Debian", "Fedora"]), "Debian")
        self.assertEqual(autoconfig.pick_distro([]), "")

    def test_gmx_share_from_windows_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            (root / "bin").mkdir()
            (root / "bin" / "gmx.exe").write_bytes(b"")
            (root / "share" / "gromacs" / "top").mkdir(parents=True)
            got = autoconfig.gmx_share_from(str(root / "bin" / "gmx.exe"))
            self.assertTrue(got.endswith("/share/gromacs/top"))

    def test_gmx_share_empty_when_missing(self):
        self.assertEqual(autoconfig.gmx_share_from(""), "")
        self.assertEqual(autoconfig.gmx_share_from("/no/such/gmx.exe"), "")

    def test_merge_only_fills_empty(self):
        cfg = {"gmx": sys.executable, "wsl_home": "", "amber_home": ""}
        det = {"gmx": "/new/gmx.exe", "wsl_home": "/home/u",
               "amber_home": "/opt/amber"}
        out = autoconfig.merge_detected(cfg, det)
        self.assertEqual(out["gmx"], sys.executable)
        self.assertEqual(out["wsl_home"], "/home/u")
        self.assertEqual(out["amber_home"], "/opt/amber")

    def test_merge_overwrite(self):
        cfg = {"gmx": "/old/gmx.exe"}
        out = autoconfig.merge_detected(cfg, {"gmx": "/new/gmx.exe"},
                                        overwrite=True)
        self.assertEqual(out["gmx"], "/new/gmx.exe")

    def test_default_config_has_no_personal_paths(self):
        for key in ("gmx", "wsl_distro", "wsl_home", "amber_home", "mmpbsa_env"):
            self.assertEqual(mdpipe.DEFAULT_CONFIG[key], "",
                             f"{key} 默认值应为空，交给自动探测")


# ---------------------------------------------------------------- 杂项
class TestMisc(unittest.TestCase):
    def test_ps_escape(self):
        self.assertEqual(notify._ps_escape("it's\nhere"), "it''s here")

    def test_size_format(self):
        self.assertEqual(mdpipe._size_h(0), "0 B")
        self.assertEqual(mdpipe._size_h(2048), "2.00 KB")


if __name__ == "__main__":
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    unittest.main(verbosity=2)
