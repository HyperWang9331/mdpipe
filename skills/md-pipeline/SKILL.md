---
name: md-pipeline
description: 用 mdpipe 命令行端到端操作 GROMACS 蛋白–配体分子动力学流水线：环境体检、投递与排队、监控、RMSD/RMSF/Rg/SASA/氢键/接触/聚类分析、MM-PBSA 结合自由能、PCA/FEL、出图与归档；遇到复杂问题会先用 mdpipe doc 检索随附的操作手册与安装清单再动手。只要用户提到分子动力学、MD、GROMACS、mdpipe、跑轨迹、投递作业、FABP5、阿魏酸–Keap1、pipeline/jobs、结合自由能或 MM-PBSA，就用这个 skill——即使他们只说“跑个模拟”“看看这批结果”。
---

# ClawsGO MD 流水线（mdpipe）

这套软件把「蛋白 + 小分子 → 可分析轨迹」的整条链路封装成一条命令：结构检查与补全、
配体参数化（AmberTools）、能量最小化、NVT、NPT、生产 MD、分析、出图、MM-PBSA、归档。
你的任务是替用户正确地驱动它，并在每一步给出可核对的结论。

## 0. 先确认入口和配置

两种形态，命令完全一样：

| 形态 | 调用方式 |
|---|---|
| 安装版 | `mdpipe …`（或 `%LOCALAPPDATA%\ClawsGO\MDpipeline\mdpipe.bat …`） |
| 源码版 | `python scripts/mdpipe.py …` |

在仓库里工作时优先用源码版，避免改了代码却跑到安装副本上。

所有机器相关的路径只放在 `pipeline/config.json`：`gmx`、`wsl_distro`、`amber_home`、
`mmpbsa_env`、`gmx_share`、`wsl_home`、`ligplot_home`、`gpu_id`、`use_gpu`、`ntomp`。
其余是方法学默认值（力场、水模型、盐浓度、盒子边距、NVT/NPT 时长）。改配置前先读一遍，
不要凭猜测改路径。

## 1. 铁律（每条都有代价）

1. **动手前先看现场**：`mdpipe doctor` + `mdpipe status`。有作业在 `running` 时不要再起
   新的 `mdrun`——单卡机器上第二个进程会把两个都拖慢，甚至 CUDA OOM。
2. **分析轨迹禁止 `-pbc nojump`**。用 `md_center.xtc`（全体系居中）或 `md_prot.xtc`
   （Protein_LIG）。首帧镜像错位会一路带下去，RMSD/RMSF 全假，而 Rg/SASA 看着正常，
   极难发现。
3. **`pipeline/jobs` 是数据，不是缓存**。清理用 `mdpipe archive`；不要 `rm -rf`，
   也不要用 `remove` 代替删除（它默认只移出队列）。
4. **分析代码在作业启动时快照**。改过 `scripts/` 里的分析逻辑后，必须手动重跑
   `mdpipe analysis`，否则旧作业仍用旧代码出图。
5. **长作业先确认**：`--ns` 大于 20、或一次投递超过 4 个作业时，先说明预计耗时和磁盘
   占用，得到用户同意再开跑。
6. **不覆盖用户已有结果**：`run --force` 会重算已完成作业，用之前必须问。

## 2. 命令地图

| 命令 | 用途 | 常用参数 |
|---|---|---|
| `doctor` | 只读环境体检 | `--plan` 只看安装清单；`--json` 机器可读；`--install-pip` 代装必需 Python 包 |
| `add` | 登记作业（一蛋白多配体、多蛋白多配体） | `--protein/--protein-dir`、`--ligand/--ligand-dir`、`--name`、`--replicates`、`--ns`、`--forcefield`、`--ligand-itp/-gro/-posre` |
| `complex` | 从「蛋白+配体」复合物 PDB 直接建作业 | `--lig-resname`、`--charge`、`--chain` |
| `scan` | 扫 `pipeline/inbox/` 自动登记 | 无参数 |
| `run` | 执行队列 | `--job`、`--ns`、`--workers`、`--gpu-jobs`、`--prep-only`、`--only-failed`、`--retry N`、`--force` |
| `status` | 看每个作业的状态/准备步数/生产 ns/自检 | 无参数 |
| `selfcheck` | 检查 xvg 数值是否合理 | `--job`、`--all` |
| `analysis` | RMSD、RMSF、Rg、SASA、氢键 + 图 | `--job`、`--all`、`--replicates`、`--force` |
| `figures` | 单体系图 / 多体系对比图 | `--job`、`--all`、`--compare`、`--labels` |
| `interactions` | 接触、氢键、盐桥、2D 示意图、LigPlot+ | `--cutoff`（nm，默认 0.45）、`--dt`（ns）、`--no-ligplot` |
| `mmpbsa` | MM-PBSA 结合自由能（WSL + gmx_MMPBSA） | `--frames`（默认 10）、`--last-ns`（默认 20）、`--keep-center` |
| `fel` | PCA / 自由能地形 | `--compare`、`--force` |
| `viz` | PyMOL 快照 + B-factor | `--frames`（ns）、`--no-pymol` |
| `archive` | 归档轨迹、校验副本 | `--copy-to`、`--analysis-only`、`--delete`（不加只预演） |
| `provenance` | 记录/回填 gmx、AmberTools、mdp 的哈希 | `--backfill` |
| `report` | 汇总到 `pipeline/summary.csv` | 无参数 |
| `remove` | 移出队列 | `--keep` 保留目录 |
| `doc` | 查操作手册与安装清单 | `--list`、`--grep 关键词`、`--section 20.10` |
| `skill` | 把本技能装到各 AI Harness | `--list`、`--target detected`、`--print` |
| `gui` | 图形界面（人用；自动化不要走 GUI） | 无参数 |

## 3. 标准流程

### 3.1 新机器接手

```bash
mdpipe doctor           # 逐项检测，标 √ / ! / ×
mdpipe doctor --plan    # 缺什么、哪些能自动装、哪些必须手动
```

「必需」项出现 × 时先补环境，不要硬跑。缺 GPU 就把 `config.json` 的 `use_gpu` 设为
`false`（纯 CPU 约慢 30 倍），不要假装 GPU 可用。

### 3.2 投递作业

一个蛋白配多个配体：

```bash
mdpipe add --protein receptor.pdb --ligand lig_1.sdf lig_2.sdf \
           --name fabp5 --replicates 2 --ns 100
```

已有复合物 PDB：

```bash
mdpipe complex "complex.pdb" --lig-resname LIG --replicates 2 --ns 100
```

批量目录：

```bash
mdpipe add --protein-dir pdb/ --ligand-dir ligs/
mdpipe scan                     # 或把文件丢进 pipeline/inbox/ 再 scan
```

配体参数化默认走 AmberTools（WSL）。如果用户已经自备参数，加
`--ligand-itp x.itp --ligand-gro x.gro --ligand-posre x.itp`，这条路径不需要 AmberTools。
GROMOS 54A7 + ATB 的配体用 `--forcefield gromos54a7 --ligand-atb …`。

### 3.3 执行

```bash
mdpipe run                       # 队列里所有未完成作业
mdpipe run --job fabp5 --ns 100
mdpipe run --prep-only           # 只做结构准备，不占 GPU
mdpipe run --retry 3 --retry-wait 120   # 笔记本/断电/驱动重置场景
```

资源由两个参数控制：`--workers`（总并发）和 `--gpu-jobs`（同时占 GPU 的 `mdrun` 数）。
默认 2 / 2；单卡且要并行做别的分析时降到 1。

### 3.4 监控

```bash
mdpipe status
mdpipe selfcheck --all
```

状态取值：`running`、`prepared`（准备完成未生产）、`done`、`stopped`（手动停）、
`failed`。`status` 还会显示开工前结构体检结论与 EM 最大力。

### 3.5 完成后分析

```bash
mdpipe selfcheck --all                     # 先确认轨迹本身没问题
mdpipe analysis --job fabp5                # 每个重复各一套曲线和图
mdpipe figures --all --compare             # 多体系叠一张对比图
mdpipe interactions --job fabp5            # 接触/氢键/盐桥 + 2D
mdpipe mmpbsa --job fabp5 --last-ns 20 --frames 10
mdpipe fel --job fabp5 --compare
mdpipe report                              # 汇总表
```

产物位置：

| 内容 | 路径 |
|---|---|
| 分析曲线与数据 | `<job>/07_analysis_r*/` |
| 单体系图 | `<job>/08_figures/` |
| MM-PBSA | `<job>/08_mmpbsa_r*/`，对比图在 `pipeline/figures/` |
| 汇总表 | `pipeline/summary.csv` |

### 3.6 归档

```bash
mdpipe archive --all --copy-to "E:\分子动力学模拟"   # 只校验副本，不复制
mdpipe archive --all --analysis-only                # 预演，不加 --delete 不删
```

## 4. 输入文件要求

蛋白：

- 单链，或用 `--chain` 明确指定；有 altloc 时用 `--force-altloc A`。
- 去水、去无关 HETATM；辅因子/金属用 `--keep-het` 显式保留。
- 缺失残基或原子先补全（`complete_protein` / 开工前体检会报）。

配体：

- SDF / MOL / PDB，必须 3D、加氢，净电荷明确（`--charge`），残基名 ≤ 3 字符。

复合物：

- 蛋白与配体在同一个 PDB，配体残基名用 `--lig-resname` 指定。

## 5. 判读结果

- EM 最大力 < 1000 kJ/mol/nm 视为正常；明显偏高说明初始结构有冲突。
- NVT 温度 300 ± 5 K；NPT 密度 1000 ± 20 kg/m³（TIP3P）。
- 生产轨迹长度应达到设定的 `prod_ns`；不足先查 `selfcheck`。
- MM-PBSA 看 ΔG 与 SEM；两个重复相差 > 5 kcal/mol 时要在结论里注明样本不足，
  不能只报更负的那个。
- 相互作用看 `occupancy`（全程占比）和 `occupancy_bound`（结合态内占比），两者含义不同。

## 6. 常见故障

| 现象 | 原因 | 处置 |
|---|---|---|
| `gmx` 找不到 / 版本报错 | `config.json` 的路径不对 | `mdpipe doctor` 确认，改 `gmx` |
| CUDA OOM、驱动重置、作业中断 | 并发太高 | `--gpu-jobs 1`，加 `--retry 3` |
| `grompp` 报 atomtype 重复 | GROMOS/ATB 的 itp 合并问题 | 见操作手册 20 章 GROMOS 线路 |
| MM-PBSA 报 BOND / 能量溢出 | 轨迹未居中或帧数异常 | 确认用 `md_center.xtc`，重跑 mmpbsa |
| 分析图空白或曲线异常 | 轨迹未完成/损坏 | `selfcheck`，必要时 `run --only-failed` |
| 作业卡在 prepared | 只跑了 `--prep-only` | 直接 `mdpipe run --job <name>` |

## 7. 复杂问题：先查手册和清单，不要猜

这份 skill 只是最短路径。完整依据在四份文档里，先用命令定位，再只读需要的那一节：

```bash
mdpipe doc --list                              # 四份文档的真实路径
mdpipe doc manual --section 16                 # 按章节号取整节
mdpipe doc manual --grep "GROMOS" --context 2  # 按关键词检索，带行号
mdpipe doc checklist --grep "AmberTools"
```

安装版的文档在 `%LOCALAPPDATA%\ClawsGO\MDpipeline\docs\`；源码版在仓库根目录、
`installer\` 和 `pipeline\`。`mdpipe doc` 找不到时，按这两个位置直接读文件。

问题与章节的对应关系：

| 问题 | 查哪里 |
|---|---|
| 输入文件不合规、从对接结果里挑文件 | manual §1.4–1.7 |
| altLoc 冲突、缺失残基、补 loop | manual §3 |
| 配体提取、加氢、净电荷 | manual §4 |
| GAFF2 + AM1-BCC 参数化 | manual §5 |
| 拓扑组装、坐标合并 | manual §6 |
| 溶剂化、加离子、索引组 | manual §7–8 |
| EM / NVT / NPT / 生产 的判据 | manual §9–12 |
| PBC、RMSD/RMSF/Rg/SASA、氢键语法、MM-PBSA 原理 | manual §13 |
| 参数速查（cutoff、温度、密度、约束） | manual §14 |
| 性能与时间预算、重复数怎么定 | manual §15 |
| 换蛋白/配体/力场、膜蛋白、糖链、共价配体 | manual §16 |
| 报错信息对照 | manual §17 |
| 批量流水线、GUI、归档、doctor | manual §20 |
| 空机器装什么、哪些能自动装 | checklist 全文；配合 `mdpipe doctor --plan` |
| 安装包行为、卸载保留什么 | install 文档；checklist §6 |

使用规则：

- 命中了上表，就先把对应节读出来再动手。力场、约束、分析参数这类改动，凭记忆改的代价
  是重跑几十小时。
- 只取需要的节。手册约 95 KB，用 `--section` / `--grep`，不要整份读进上下文。
- 回答里写清章节号（例如「依据手册 §16 换力场」），用户能自己复核。
- skill 与手册冲突时以手册为准：skill 是速查，手册是依据。
- 用户问“这台机器还缺什么”“能不能离线”，先看 checklist 和 `doctor --plan`，那里把
  「能自动装 / 半自动 / 只能手动」分得很清楚。

## 8. 向用户汇报

每次操作后给一段可直接核对的结论，至少包含：作业名、状态、生产时长、自检结论、
关键指标（EM 最大力、温度、密度、ΔG 等）、产物路径。不要只把日志贴出来。

## 9. 不要做

- 不要在没有确认的情况下删作业、改机器相关路径、或重算已完成作业。
- 不要在分析里加 `-pbc nojump`。
- 不要用 GUI 点击代替 CLI 自动化；GUI 是给人用的，命令行才是可复现的。
- 不要把绝对路径写进脚本或命令历史；需要新路径就写进 `config.json`。
