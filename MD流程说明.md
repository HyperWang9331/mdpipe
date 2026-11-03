# 本机分子动力学（蛋白 + 小分子）运行方案

> 环境实测日期：2026-10-29 ｜ 工作机：Windows 11 + RTX 5060 Laptop (8 GB) ｜ 所有性能数字均为本机实测，非估算

---

## 0. 一句话结论

本机 **GROMACS 2026.3 原生 CUDA 版已就绪，GPU 加速可用**，直接跑蛋白+小分子的动力学没有问题；唯一缺的一环是**小分子的力场参数化工具**（Windows 上没有 AmberTools），需要在下面三条路线里选一条。给我蛋白结构 + 配体结构，我可以从建拓扑一路做到分析和出图。

---

## 1. 本机实测环境

| 项目 | 实测结果 | 说明 |
|---|---|---|
| GROMACS | **2026.3**（`<GROMACS目录>\bin\gmx.exe`），mixed precision，thread-MPI | 已加入 PATH，直接 `gmx` 调用 |
| CUDA | 编译目标 **sm_120**，驱动 13.40 | 原生支持 RTX 5060 (Blackwell)，无需 JIT 回退 |
| GPU | RTX 5060 Laptop，8 GB 显存 | 显存是体系规模的上限因素 |
| CPU | 24 逻辑核，AVX2_256 | |
| 磁盘 | C: 剩 36 GB / D: 剩 36 GB（均 92-93% 占用） | **轨迹文件的主要约束** |
| 内置力场 | amber03/14sb/19sb/99sb-ildn、charmm27、oplsaa、gromos43a1–54a7 | 无 CHARMM36m、无 Martini |
| 配体工具 | **无 antechamber / acpype / tleap**（win-64 无 AmberTools 包） | 见第 4 节三条路线 |
| 辅助工具 | OpenBabel 3.1.1、RDKit (py3.14)、PyMOL | 格式转换 / 构象生成 / 看图 |
| WSL | Ubuntu-24.04 已装且可用，Python 3.12，**未装 conda/AmberTools** | 路线 A 的落点 |
| 网络 | PyPI 直连不通，清华镜像可用；GitHub、CHARMM-GUI、LigParGen 可达；ambermd.org / mackerell.umaryland.edu 不通 | 装包一律走清华镜像 |

### 实测吞吐量（纯水体系，PME + LINCS，dt = 2 fs）

| 体系规模 | 原子数 | GPU（`-ntomp 4`） | 纯 CPU（12 线程） |
|---|---|---|---|
| 4.2 nm 立方 | 7,362 | **576 ns/day** | — |
| 7.0 nm 立方 | 33,675 | **346 ns/day** | 39.7 ns/day |
| 8.4 nm 立方 | 57,801 | **230 ns/day** | — |

由此推算一个典型的 5 万原子蛋白–配体复合物：

| 目标时长 | GPU 墙钟时间 |
|---|---|
| 100 ns | ≈ 11 小时 |
| 500 ns | ≈ 2.2 天 |
| 1 µs | ≈ 4.3 天 |

**调优结论（实测）**：`-ntomp 4` 最快；给到 8 或 16 线程不再有增益（8 GB 显存 + 内存带宽瓶颈，多出来的核心在空转）。`-bonded gpu` 可用：19,430 原子体系上 5 万步实测 455–480 ns/day，比不加时的 432 ns/day 快 5%~11%。

推荐命令：

```bash
gmx mdrun -deffnm production -ntmpi 1 -ntomp 4 -gpu_id 0 -nb gpu -pme gpu -bonded gpu
```

---

## 2. 整体流程（八步）

```
① 结构预处理  →  ② 蛋白拓扑  →  ③ 配体拓扑  →  ④ 组装体系
   ↓
⑤ 溶剂化+离子  →  ⑥ 能量最小化  →  ⑦ NVT/NPT 平衡  →  ⑧ 生产 MD → 分析
```

### ① 结构预处理

```bash
# 检查 PDB：缺失残基、交替构象、非标准残基、结晶水/金属离子
python -c "import rdkit; print(rdkit.__version__)"   # 需要时用 RDKit 检查配体
```

要确认的事：缺失原子是否补全（pdb2gmx 会补氢，但不补主链缺口）、His 质子化态、链 ID、是否保留晶体水/金属离子/辅因子、是否要做突变。

### ② 蛋白拓扑

**关键**：pdb2gmx 不认识小分子残基，必须先把配体从 PDB 里摘出去。

```bash
# 仅蛋白部分 → 拓扑 + 坐标 + 位置约束
echo "1 1 1" | gmx pdb2gmx -f protein.pdb -o protein.gro -p topol.top \
    -i posre.itp -ff amber14sb -water tip3p -ignh -ss
#        ↑ His 质子化态: 0=HID 1=HIE 2=HIP, 每个 His 一个数字
```

- `-ignh`：丢弃晶体里的氢，按力场重新加（晶体氢的位置和质子化态不可靠）
- His 状态是最容易出错、影响又最大的一个选择——活性位点的 His 猜错，氢键网络整个变样。建议先用 PyMOL 检查附近氢键受体再定
- 输出 `topol.top` / `protein.gro` / `posre.itp`

### ③ 配体拓扑（唯一需要你选路线的环节）

三条路线见第 4 节，产物都是一个 `lig.itp`（含 `[ atoms ]` 电荷、`[ bonds ]`、`[ dihedrals ]`）和可选的 `posre_lig.itp`。

### ④ 组装体系

把蛋白坐标和配体坐标合成一个 `.gro`，再把配体分子类型挂进 `topol.top`：

```python
# 合并坐标：读 protein.gro + lig.gro，原子顺序 = 蛋白在前、配体在后
import sys
def read_gro(p):
    L = open(p).read().splitlines()
    n = int(L[1]); return L[2:2+n], L[2+n].split()

patoms, pbox = read_gro("protein.gro")
latoms, _    = read_gro("lig.gro")
out = ["complex", f"{len(patoms)+len(latoms):5d}"]
out += patoms + latoms
out.append("   ".join(f"{float(x):9.5f}" for x in pbox[:3]))
open("complex.gro","w").write("\n".join(out)+"\n")
```

`topol.top` 的结构必须是这个顺序：

```
#include "amber14sb.ff/forcefield.itp"     ; 必须最先
#include "lig.itp"                          ; 配体分子类型（自定义原子类型在此定义）
#ifdef POSRES
#include "posre.itp"                       ; 蛋白位置约束
#include "posre_lig.itp"                   ; 配体位置约束（若有）
#endif
#include "amber14sb.ff/tip3p.itp"          ; 与力场配套的水模型
#include "amber14sb.ff/ions_tip3p.itp"     ; 离子（注意是 ions_tip3p.itp，不是老的 ions.itp）

[ system ]
Protein–ligand complex

[ molecules ]
Protein_chain_A    1
LIG                1
SOL             9273      ; 后续 solvate/genion 自动更新
NA                 8
CL                 0
```

### ⑤ 溶剂化 + 加离子

```bash
gmx editconf -f complex.gro -o box.gro -bt dodecahedron -d 1.0 -c
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

# genion 需要一个 .tpr
gmx grompp -f mdp_templates/em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 0
echo "SOL" | gmx genion -s ions.tpr -o system.gro -p topol.top \
     -pname NA -nname CL -neutral -conc 0.15
```

`-d 1.0` 表示蛋白到盒子边缘 1.0 nm（对周期镜像即最小间距 2.0 nm，配合 1.0 nm 截断足够）。菱形十二面体比立方盒子省约 29% 的水分子，是球状蛋白的默认选择。

### ⑤' 建立索引组（控温分组需要）

`nvt/npt/production.mdp` 里用的是 `Protein_LIG` 和 `Water_and_ions` 两个耦合组。`Protein_LIG` 要自己建，`Water_and_ions` 是 GROMACS 2023 以后的默认组（本次第 20 组），不用手工拼：

```bash
printf "1 | 13\nq\n" | gmx make_ndx -f system.gro -o index.ndx
# 1 = Protein, 13 = LIG（编号以 gmx make_ndx 屏幕输出为准）
# 并集会被自动命名为 Protein_LIG，不要再写 name，否则组号错位
```

（`16 | 19` 拼出来的组名是 `Water_Ion`，不是 `Water_and_ions`；既然默认组已有现成的，就别多这一步。）

之后每一阶段的 `grompp` 都加 `-n index.ndx`。耦合组必须覆盖所有原子，漏掉配体或离子会让 grompp 报 "atom not covered by a temperature coupling group"。

### ⑥⑦⑧ 最小化 → 平衡 → 生产

```bash
# ⑥ 能量最小化（EM 是单组控温，不需要 -n）
gmx grompp -f mdp_templates/em.mdp -c system.gro -p topol.top -o em.tpr -maxwarn 0
gmx mdrun -deffnm em -ntmpi 1 -ntomp 8
grep "Maximum force" em.log | tail -1     # 必须 < 1000 kJ/mol/nm

# ⑦ NVT（-r 提供位置约束参考坐标）
gmx grompp -f mdp_templates/nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 0
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 4 -gpu_id 0 -nb gpu -pme gpu -bonded gpu

# ⑦ NPT（-t 继承 NVT 速度）
gmx grompp -f mdp_templates/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 0
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 4 -gpu_id 0 -nb gpu -pme gpu -bonded gpu

# ⑧ 生产（无 -DPOSRES，约束解除）
gmx grompp -f mdp_templates/production.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr -maxwarn 0
gmx mdrun -deffnm md -ntmpi 1 -ntomp 4 -gpu_id 0 -nb gpu -pme gpu -bonded gpu -maxh 23.5
```

**每个阶段都必须过的检查**：

| 阶段 | 判据 |
|---|---|
| grompp | `-maxwarn 0` 且无 WARNING；有 WARNING 一律先查清再放行 |
| EM | `Maximum force < 1000`；势能为明显负值 |
| NVT | 温度在 20–30 ps 内稳定到 300 K；无 LINCS warning |
| NPT | 密度落在所选定 `DispCorr` 对应的区间：`EnerPres` 为 1020–1040 kg/m³（实测 1030.9），`no` 为 1010–1030（实测 1019.5）；压力均值 1 bar |
| 生产前 | 用 panedr 看 T/P/密度/能量的分块均值是否平稳，不平稳就延长平衡 |

---

## 3. 我的默认模拟参数

**力场组合：AMBER ff14SB（蛋白）+ GAFF2（配体）+ TIP3P（水）+ PME + 1.0 nm 截断**。理由：GAFF2 与 ff14SB 同为 AMBER 家族，1-4 缩放一致（fudgeLJ 0.5 / fudgeQQ 0.8333），配体参数与蛋白参数天然兼容；这是当前蛋白–小分子模拟最成熟、文献最多的组合。

模板已放在 `mdp_templates/`：`em.mdp`、`nvt.mdp`、`npt.mdp`、`production.mdp`。

| 参数 | 取值 | 为什么 |
|---|---|---|
| integrator | `md`（leap-frog），dt = 2 fs | X–H 键受约束后 2 fs 稳定；4 fs 需要 all-bonds + 更高 LINCS 阶数，收益有限 |
| constraints | `h-bonds` + LINCS (order 4, iter 1) | 冻结最快的振动，允许 2 fs 步长 |
| 静电 | PME，rcoulomb = 1.0 nm，fourierspacing 0.16 nm，pme-order 4 | AMBER 惯例 1.0 nm；生物体系必须 PME，纯截断会产生假象 |
| 范德华 | rvdw = 1.0 nm，`vdw-modifier = Potential-shift`，`DispCorr = EnerPres` | AMBER 力场用 shift 而非 force-switch；**漏掉 DispCorr 会让 TIP3P 密度偏低约 1%** |
| 盒子 | 菱形十二面体，蛋白距边缘 1.0 nm | 体积比立方盒子小约 29%，电荷体系建议 1.2 nm |
| 盐 | 0.15 M NaCl（`-neutral -conc 0.15`） | 生理条件；纯理论对照可用 0 |
| NVT | 500 ps，V-rescale (τ=0.1 ps)，`-DPOSRES` | 正则系综正确；位置约束让溶剂先弛豫、蛋白不动 |
| NPT | 1 ns，C-rescale (τ-p = 2 ps, 1.013 bar, κ=4.5e-5) | C-rescale 收敛快且系综正确，冷启动安全；Parrinello-Rahman 从非平衡密度启动会剧烈振荡 |
| 生产 | V-rescale；如需严格 NPT 涨落再切 Parrinello-Rahman (τ-p ≥ 2 ps) | V-rescale 对结构/热力学平均足够；要算热容或速度分布则换 Nose-Hoover |
| 非键参数 | 四个阶段完全一致 | 阶段间截断/PME 设置不一致是最隐蔽的错误来源 |
| 轨迹输出 | 每 20 ps 一帧 xtc（生产） | 58k 原子 ≈ 0.1 MB/帧 → 1 µs 约 5 GB；磁盘只剩 36 GB，可按需放宽到 50 ps |

**如果研究目标是配体–蛋白结合自由能**，还需要在平衡后加一段自由能计算（GROMACS 的 lambda 耦合 + BAR/MBAR），或者用 MM-PBSA 做粗略排序（需要额外安装 gmx_MMPBSA）。这部分要不要做、怎么做，取决于你的科学问题，可以单独讨论。

---

## 4. 小分子参数化：三条路线

### 路线 A（推荐，本地全自动）：WSL 里装 AmberTools，用 antechamber + acpype 生成 GAFF2

```bash
# 一次性配置（WSL Ubuntu-24.04）
wsl -d Ubuntu-24.04
curl -L https://mirrors.tuna.tsinghua.edu.cn/github-release/conda-forge/miniforge/... # 或直接装 Miniforge
conda create -n amber -c conda-forge ambertools -y     # linux-64 有 ambertools 26.0，约 3 GB
conda activate amber

# 每个配体
acpype -i lig.mol2 -b lig -c bcc -a gaff2 -n <净电荷>     # 产物: lig.itp, posre_lig.itp, lig.gro
```

- 优点：可脚本化、可重复、批量处理多个配体；AM1-BCC 电荷（BCC）质量可靠
- 代价：一次性约 3 GB 磁盘 + 半小时配置；WSL 和 Windows 的 GROMACS 分工明确（WSL 只做拓扑，MD 在 Windows GPU 上跑）
- 前提：conda-forge 的 ambertools **只有 linux-64/osx 版本，没有 win-64**，所以必须走 WSL

### 路线 B（零安装）：LigParGen 网页 → OPLS-AA 参数

上传配体 mol2/PDB，网站直接返回 GROMACS 格式的 `.itp` + `.gro`；蛋白侧改用 `gmx pdb2gmx -ff oplsaa -water tip4p`。

- 优点：不用装任何东西，几分钟出结果
- 代价：OPLS-AA 蛋白力场比 ff14SB 老一些；需要上传结构到第三方网站；一个配体一次，不适合批量

### 路线 C（最稳，质量最高）：CHARMM-GUI Ligand Reader & Modeler

上传蛋白 + 配体，选 CHARMM36m + CGenFF，网站直接给出**完整可跑的 GROMACS 输入包**（拓扑、mdp、索引组、运行脚本），下载后本机 `mdrun` 即可。

- 优点：拓扑问题最少，膜蛋白/糖基化/共价配体都支持；CGenFF 对药物分子覆盖好
- 代价：需要注册账号、网页交互，我无法代你完成上传；本机没有 CHARMM36m 力场文件（需从镜像获取），但 CHARMM-GUI 的包里自带

**我的建议**：

- 常规可溶蛋白 + 药物样小分子 → **A**。一次配置（约 30 分钟、3 GB 磁盘），之后完全由我在本机自动完成，不需要你上传结构到任何网站，也能批量处理几十个配体做构效关系。
- 体系复杂时（膜蛋白、糖基化、共价配体、金属配位配体、非标准残基多）→ **C**。CHARMM-GUI 的拓扑模板质量最高，代价是你要在网页上操作一遍。
- 只想快速试一下、不想装任何东西 → **B**，但 OPLS-AA 蛋白力场偏老，正式研究不建议。

---

## 5. 需要你提供的文件与信息

### 必须提供

| # | 内容 | 格式 | 说明 |
|---|---|---|---|
| 1 | **蛋白结构** | `.pdb`（或 `.cif`） | X-ray / cryo-EM / AlphaFold 均可。若是复合物晶体结构，最好给**原始 PDB**，配体在里面我更容易对齐 |
| 2 | **小分子结构** | `.sdf` / `.mol2` / SMILES | 带三维坐标最好；只给 SMILES 我用 RDKit 生成构象。**务必写明净电荷** |
| 3 | **研究问题** | 文字 | 想看稳定性？结合模式？构象变化？比较几个配体？——决定模拟时长和分析方案 |

### 强烈建议一并说明

| # | 内容 | 为什么重要 |
|---|---|---|
| 4 | **配体在结构中的位置** | 晶体结构里已有配体最好；若需对接，告诉我结合口袋残基 |
| 5 | **质子化状态 / pH** | His 选 HID/HIE/HIP，Asp/Glu 是否带电，配体有无可电离基团——猜错会改变整个氢键网络 |
| 6 | **保留哪些非蛋白组分** | 结晶水、金属离子（Zn²⁺/Mg²⁺）、辅因子、结构域间的二硫键 |
| 7 | **缺失残基/loop 处理** | 补全（Modeller/Homology）还是截断保留缺口 |
| 8 | **突变体信息** | 若结构与实际研究对象不同 |
| 9 | **力场偏好** | 默认 ff14SB + GAFF2 + TIP3P；若有文献要求 CHARMM36m/CGenFF 或 OPLS-AA 请说明 |
| 10 | **盐浓度 / 温度** | 默认 0.15 M NaCl、300 K |
| 11 | **模拟时长与重复次数** | 默认 500 ns × 1 条；发表级结论建议 3 个独立种子各 ≥ 500 ns |
| 12 | **磁盘预算** | 目前只剩 36 GB，1 µs 轨迹约 5 GB（20 ps/帧）。空间紧张就降低存帧频率 |

### 需要你决定的一件事

**配体参数化路线选 A / B / C 哪一条**（第 4 节）。这决定我要不要先在 WSL 里装 AmberTools（约 30 分钟、3 GB 磁盘）。

---

## 6. 我会交付什么

```
project-20261029-1/
├── 00_structure/      protein.pdb, lig.sdf（你提供）+ 预处理记录
├── 01_topology/       topol.top, lig.itp, posre*.itp, complex.gro, system.gro
├── 02_mdp/            em.mdp, nvt.mdp, npt.mdp, production.mdp（已提供模板）
├── 03_em/  04_nvt/  05_npt/  06_production/
├── 07_analysis/       xvg + 图
└── REPORT.md          参数记录、grompp 日志、.tpr 的 SHA-256、分析与结论
```

标准分析会包含：主链/配体 RMSD、每残基 RMSF、回旋半径、SASA、蛋白–配体氢键占有率、接触频率（MDAnalysis）、关键距离随时间演化；需要的话再加 PCA / 自由能景观（FEL）。

---

## 7. 已知限制与风险

1. **无 AmberTools（Windows）**：GAFF2 参数化必须走 WSL 或网页，这是唯一必须解决的卡点。
2. **磁盘 36 GB**：1 µs 生产轨迹 + 中间文件约 8–10 GB，可以跑，但建议 1 µs 以内单次，或降低存帧频率、及时清理中间 `.trr`。
3. **8 GB 显存**：5 万原子以内很舒服；20 万原子以上（膜蛋白、大复合物）吞吐量会明显下降，可能需要考虑更大体系或云端资源。
4. **`-bonded gpu` 可用但非必需**：19k 原子体系实测 +5%~11%（455–480 vs 432 ns/day），推荐加上；若某个体系加了报错，去掉即可。
5. **力场库不全**：无 CHARMM36m、无 Martini。要跑膜体系或粗粒化，需先补力场文件（可从 GitHub 镜像取，官网不通）。
6. **网络限制**：PyPI 直连不通，pip 装包请加 `-i https://pypi.tuna.tsinghua.edu.cn/simple`；ambermd.org 与 mackerell.umaryland.edu 不可达。
7. **Windows 路径**：工程目录不要含中文和空格（当前目录符合要求）。GROMACS 在 Windows 上对超过 260 字符的长路径偶尔出错，目录层级保持浅。
8. **单机单卡**：没有 MPI 集群，多副本并行会互相抢 GPU，建议串行跑。

---

## 8. 下一步

1. 你把上面的结构文件发给我，并告诉我研究问题与配体参数化路线（A/B/C）
2. 我先做结构检查和预处理，把不确定的地方（缺失残基、His 状态、配体位置）列出来跟你确认
3. 建拓扑 → 最小化 → 短平衡（1 ns 级别）**做一次快速验证**，确认能量、密度、RMSD 正常后再提交长程生产
4. 生产跑完做标准分析并出报告
