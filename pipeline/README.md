# pipeline/ —— 自动流水线工作区

由 `scripts/mdpipe.py` 使用，对应《分子动力学操作手册》的十一步。手工流程仍在手册里；这里是"投文件、排队、出结果"的快速通道。

## 要投哪两个文件

不管对接用的是什么软件，最后都只要两份：**只有蛋白的 PDB** + **只有配体的位姿**（`.pdbqt`/`.mol2`/`.sdf`）。

- Vina / AutoDock：蛋白传对接时用的**原始受体 PDB**（`receptor.pdbqt` 不行），配体传 `output.pdbqt`。
- 只拿到"蛋白+配体"复合物 PDB（CB-Dock 之类）：**有专门的复合物入口，自动拆分**。

  ```bash
  python scripts/mdpipe.py complex complex.pdb --ns 100            # 命令行，一个文件
  python scripts/mdpipe.py complex a.pdb b.pdb c.pdb --ns 100      # 一次多个
  # 或者：inbox/Gal_1/ 里只放复合物，然后 scan
  # 或者：界面里「复合物 PDB」那一行，选文件（可多选）→ 添加复合物
  ```

  想先看看拆出什么，也可以单独跑 `python scripts/split_complex.py --in complex.pdb --out-dir split`。
  拆出来的两份存在 `pipeline/splits/<作业名>/`，原始复合物路径记在 `job.json` 的 `source_complex`。

详细对照表见《分子动力学操作手册》§1.5。

## 三步开始

```bash
# 1) 把文件丢进 inbox/<作业名>/：一个蛋白 .pdb + 一个（或多个）配体 .pdbqt
#    一个子目录放多个配体时，每个配体各建一个作业
#    只放一个"蛋白+配体"复合物 .pdb 时，scan 会自动拆开
#    需要非默认净电荷时，在同目录写一个 charge.txt（内容只有数字，如 -1）
mkdir -p inbox/Gal_1
cp FABP5.pdb output.pdbqt inbox/Gal_1/

# 2) 扫描并登记
python scripts/mdpipe.py scan

# 3) 跑（workers=同时处理的作业数，gpu-jobs=同时占用 GPU 的 mdrun 数；默认都是 2）
python scripts/mdpipe.py run
```

也可以用一条命令直接登记：

```bash
python scripts/mdpipe.py add --protein FABP5.pdb --ligand output.pdbqt \
    --name Gal_1 --charge 0 --ns 200
```

## 开工前的自动体检

投进去的蛋白不一定是完整的。数据库和 PyMOL 清洗过的受体常缺整段 loop，**pdb2gmx 遇到主链缺口不报错**——它把一条链拆成两条、在断口补上带电末端，拓扑看着正常，体系从第一步就错了。

所以「protein」这一步现在先体检再动手：

| 体检结论 | 流水线怎么处理 |
|---|---|
| PASS | 照旧走 `prep_protein.py` |
| WARN（缺侧链）/ FAIL（缺 loop、断链、修饰残基） | 先跑 `prep_receptors_full.py` 补全，用 `00_input/_prep/<名>_ready.pdb` 继续，原文件不动 |
| 体检本身出错 | 记一条警告，按原流程继续 |

结论记在 `job.json` 的 `structure_check` 里，`mdpipe status` 的「体检」一列能看到，补过的会标成 `FAIL·已补`。详细原因在作业目录的 `logs/structure_check.log`，补全过程在 `logs/autoprep.log`。

三种情况不自动补：

- 投递时加了 `--keep-het`（比如要留结构锌）——补全会把杂原子当小分子删掉；
- 蛋白是你自己处理过的（`--prepared`）——只体检、只记录；
- `config.json` 里设了 `"auto_repair": false`。

受体文件里连 SEQRES 都没有时（PyMOL / AutoDock 清洗过的 `*_clean.pdb` 都是这样），补链需要全长序列：

```bash
# 从 UniProt 下载，文件名要和蛋白同名
curl -s "https://rest.uniprot.org/uniprotkb/P31749.fasta" -o seq/AKT1.fasta

python scripts/mdpipe.py add --protein AKT1_clean.pdb --ligand fa.pdbqt \
    --name FA_AKT1 --seq-dir seq/
```

`--seqres-dir` 用于有原版 PDB、只是想借它 SEQRES 记录的情况。


## 图形界面（不想敲命令就用这个）

双击桌面上的「分子动力学流水线」快捷方式，或在项目根目录执行：

```bash
python scripts/mdpipe_gui.py
```

`python scripts\...` 是按当前目录解析的，在别的目录下要先 `cd` 到项目根目录；
或者用脚本的绝对路径 `python <项目>\scripts\mdpipe.py gui`，从任何位置都能启动。

界面分三块：

- **新建作业**：第一行放只有蛋白的 PDB，「复合物 PDB」那一行一次可选多个蛋白+配体合在一起的文件（点「添加复合物」自动拆），第二行放配体位姿（可 Ctrl 多选，也可选整个文件夹）。填作业名/净电荷/生产时长/重复数，点「添加蛋白+配体」。净电荷勾「自动判定」就交给程序推。
- **作业列表**：每个作业的状态、准备步骤进度、NPT 密度、ns/day。
- **运行日志**：实时滚动，和命令行输出一模一样。

按钮：`批量配对…` 选两个文件夹做「多蛋白 × 多配体」；`开始运行` 跑队列；`停止` 让当前这一步做完后停下（已完成的步骤保留，随时可续）；`打开作业目录` 直接进资源管理器；`导出汇总` 生成 summary.csv；`出图`（Ctrl 多选则出对比图）；`清空历史`（移出已完成的作业，问你要不要连目录一起删）。下面一条进度条显示已完成步骤数和剩余时间估算。

图形界面只是 `mdpipe.py` 的壳，两者可以混用：界面里加作业、命令行跑，或者反过来都行，进度是按作业目录里的标记算的。

**菜单栏**：

- `设置 → 组件路径与参数…`：在界面里直接改 `gmx`、WSL 发行版、AmberTools / gmx_MMPBSA 环境、GMXLIB、LigPlot+ 根目录，保存即写回 `pipeline/config.json`；`恢复默认` 只是填回默认值，点保存才生效。GROMACS 路径填错会先提示再决定是否保存。
- `设置 → 运行环境自检`：等价于 `mdpipe doctor`，结果滚动到下方日志区。
- `设置 → 打开 config.json`：其余参数（力场、水模型、并发数等）直接在编辑器里改。
- `帮助 → 关于`：版本、引用的组件与文献出处。
- `帮助 → 操作手册`：打开手册。

## 常用命令

| 命令 | 作用 |
|---|---|
| `mdpipe.py doctor` | 环境自检：GROMACS / GPU / WSL / AmberTools / 出图工具 / 凭据 / 磁盘 |
| `mdpipe.py status` | 看每个作业跑到第几步（含体检、自检两列） |
| `mdpipe.py selfcheck --all` | 复检所有作业的分析结果 |
| `mdpipe.py archive --all` | 列出可重算的轨迹副本；加 `--delete` 才删 |
| `mdpipe.py viz --job Gal_1` | 复合物 RMSD / 埋藏 SASA / 质心距离 + PyMOL 快照与 B-factor 图 |
| `mdpipe.py interactions --job Gal_1` | 接触频率 / 氢键配对 / 配体聚类 / 2D 示意图（可选 LigPlot+） |
| `mdpipe.py mmpbsa --job Gal_1` | MM-PBSA 结合自由能（批量对比图在 `pipeline/figures/`） |
| `mdpipe.py report` | 汇总成 `summary.csv`（EM 最大力、NVT 温度、NPT 密度、ns/day） |
| `mdpipe.py analysis --job Gal_1` | 按当前代码重算分析（轨迹修正 + 指标 + PCA + 图），不重跑 MD |
| `mdpipe.py figures --job Gal_1` | 重出六联图（跑完会自动出，这个用来补） |
| `mdpipe.py figures --job Gal_1 Glc_1` | 多体系对比图 → `pipeline/figures/` |
| `mdpipe.py fel --job Gal_1` | PCA 自由能景观（≥10 ns 自动算，`--force` 强制） |
| `mdpipe.py fel --job Gal_1 Glc_1 --compare` | 多体系 FEL 3D 对比图 |
| `mdpipe.py run --job Gal_1` | 只跑指定作业 |
| `mdpipe.py run --job Gal_1 --prep-only` | 只做结构准备（拓扑/溶剂/离子/索引），不占 GPU |
| `mdpipe.py run --only-failed` | 只重跑失败的 |
| `mdpipe.py run --force` | 清掉完成标记，从头重跑 |
| `mdpipe.py remove Gal_1` | 从队列移除（加 `--keep` 保留目录） |

## 并发是怎么排的

- **准备阶段并行**：结构清洗、pdb2gmx、溶剂化、建索引都是 CPU 活，`--workers` 个作业同时做。
- **MD 阶段排队**：每次 `mdrun` 前要拿 GPU 名额，名额数 = `--gpu-jobs`（默认 2）。
- **默认 `--gpu-jobs 2`**：本机用 19,430 原子体系实测，两个 mdrun 同时跑（各 50,000 步）得到 227 + 234 = **461 ns/day** 合计，单个跑是 455–480 ns/day，基本持平。也就是说同时跑两个不亏，队列消化得更整齐；用户实测两个 100 ns 作业并行是效率最优点。想更省内存或单作业延迟优先，就设回 1。
- **落单作业的重复并行**：待跑作业数少于 `--workers` 时（典型情况是队列只剩最后一个体系），该作业的两个重复会同时跑——单个 mdrun 只吃 50–60% 的 GPU，两个一起能到约 540 ns/day 合计。准备步骤和 EM 仍然只做一次，重复之间只共享只读文件。`--parallel-replicates` / `--serial-replicates` 可以强制开关；某个重复已经跑了一部分时，串行通常更快。


## 一对多 / 多对多

| 场景 | 命令 |
|---|---|
| 1 蛋白 × N 配体 | `mdpipe.py add --protein P.pdb --ligand a.pdbqt b.pdbqt c.pdbqt` |
| 1 蛋白 × 一整个配体文件夹 | `mdpipe.py add --protein P.pdb --ligand-dir <配体目录>` |
| M 蛋白 × N 配体（全组合） | `mdpipe.py add --protein-dir <蛋白目录> --ligand-dir <配体目录>` |
| 只跑指定组合（命令行） | `mdpipe.py add --pair "FABP5::Glc,Gal" --pair "FABP6::Ara" --protein-dir <蛋白> --ligand-dir <配体>` |
| 只跑指定组合（表格） | `mdpipe.py add --csv pairs.csv`（列 `protein,ligand`，可选 `name,charge,ns,replicates`） |
| 界面 | 配体多选 / 选文件夹；「批量配对…」弹出全组合勾选表，点掉不要的；「导入配对表…」读 CSV |

`--pair` 用两个冒号分隔蛋白和配体（`蛋白::配体1,配体2`），因为单个冒号会和 `C:` 盘符冲突；两侧都可以写全路径或只写文件名。

作业命名：1×1 用配体名；1×N 用配体名；多个蛋白用 `蛋白名__配体名`；重名自动加 `_2`。
超过 50 个作业要加 `--yes`，防止误展开。

同一个 `output.pdbqt` 里的多个构象属于同一个分子，用 `--model 2` 选第几个，不要按多分子处理。

## 改条件：改哪儿

| 想改什么 | 怎么做 |
|---|---|
| 换蛋白 / 换配体 | 新建作业；一对多、多对多见上文 |
| 生产时长 | 添加时填「生产时长」（界面）或 `--ns 500` |
| 独立重复条数 | 填「独立重复」或 `--replicates 3`（会自动用 1001/1002/1003 三个种子） |
| 净电荷 | 填数字，或勾「自动判定」 |
| 用第几个对接构象 | 命令行 `--model 2` |
| 蛋白只保留某条链 / altLoc | `--chain A`、`--force-altloc 120:B,127:B` |
| 盐浓度、盒子边距、力场、水模型、默认并发数 | 编辑 `pipeline/config.json` |
| 开工前是否自动补全残缺蛋白 | `config.json` 的 `"auto_repair": true/false` |
| 换力场/水模型 | `config.json` 的 `"forcefield"`、`"water"`，或投递时 `--forcefield oplsaa --water spc` |
| 用自己算好的配体参数（LigParGen / ATB / CGenFF） | 加 `--ligand-itp X.itp --ligand-gro X.gro`，跳过 antechamber |
| 走 GROMOS 54A7 线路（配体从 ATB 下载） | `--forcefield gromos54a7 --ligand-atb ATB导出目录/`，自动去重 atomtype、自动换 GROMOS 的 mdp |
| 截断半径、控温控压、DispCorr 等 mdp 参数 | 编辑 `mdp_templates/*.mdp`，对所有新作业生效 |
| 分析项目 | 改 `scripts/mdpipe.py` 里的 `analysis_step` |

`config.json` 首次运行自动生成，组件路径由 `scripts/autoconfig.py` 探测后写入（安装程序也会跑一次）。改完直接生效，不用重启：

```json
{
  "gmx": "<GROMACS目录>/bin/gmx.exe",
  "wsl_distro": "Ubuntu-24.04",
  "amber_home": "/home/<用户名>/miniconda3/envs/amber",
  "mmpbsa_env": "/home/<用户名>/miniconda3/envs/mmpbsa",
  "gmx_share": "",
  "wsl_home": "/home/<用户名>",
  "ligplot_home": "",
  "gpu_id": 0, "ntomp": 4, "ntomp_em": 8,
  "forcefield": "amber14sb", "water": "tip3p",
  "salt_conc": 0.15, "box_distance": 1.0,
  "nvt_ps": 500, "npt_ps": 1000, "prod_ns": 200,
  "gen_seed_base": 1000, "gpu_jobs": 2, "workers": 2,
  "fel_min_ns": 10, "throughput_ns_per_day": 460,
  "auto_repair": true, "notify": true, "notify_per_job": false
}
```

## 环境自检与路径

```bash
python scripts/mdpipe.py doctor               # 或安装后的 mdpipe doctor
python scripts/mdpipe.py doctor --plan        # 只打印安装清单：哪些能自动装、哪些必须手动
python scripts/mdpipe.py doctor --install-pip # 代装缺失的 numpy/scipy/matplotlib（需网络）
python scripts/mdpipe.py doctor --json        # 给安装程序 / CI 用
python scripts/mdpipe.py config --detect      # 只看自动探测到的组件路径
python scripts/mdpipe.py config --detect --write  # 写入 config.json（只补空项）
```

完整的安装清单与"空电脑按顺序装"的步骤见 [安装与环境清单.md](../安装与环境清单.md)。

只读检测 Python 包、GROMACS（版本与 CUDA 支持）、GPU、WSL、AmberTools、gmx_MMPBSA、PyMOL、LigPlot+、Ghostscript、OpenBabel 和磁盘空间，缺什么给出下一步；退出码非 0 表示有必需项缺失。

**所有外部组件的位置都以 `config.json` 为准**，检测到的每一项都会注明来源：

| 字段 | 指向 | 找不到时 |
|---|---|---|
| `gmx` | `gmx.exe` 完整路径 | 再查 `PATH`；都没有就报错 |
| `wsl_distro` | WSL 发行版名 | 列出已有发行版供对照 |
| `amber_home` | AmberTools 环境（`bin/antechamber`） | 报 WARN |
| `mmpbsa_env` | gmx_MMPBSA 环境（`bin/gmx_MMPBSA`） | 报 WARN |
| `gmx_share` | WSL 里的 `share/gromacs/top`；留空则从 `gmx` 推出 | — |
| `wsl_home` | WSL 用户家目录 | — |
| `ligplot_home` | LigPlot+ 根目录；留空则查 `LIGPLUS_HOME`，再查 `<LigPlot+目录>` 等常见位置 | 跳过 LigPlot+，用 RDKit 图 |

换机器只需要改这一处，脚本里没有硬编码路径。`config.json` 首次运行自动生成，改完下次调用生效。

## 复合物可视化（viz）

六联图回答"蛋白稳不稳"，`viz` 回答"配体在不在、埋得深不深"：

```bash
python scripts/mdpipe.py viz --job Gal_1                 # 全部重复
python scripts/mdpipe.py viz --job Gal_1 --frames 0 20 40 60 80 100
python scripts/mdpipe.py viz --job Gal_1 --no-pymol      # 只算曲线，不出结构图
```

产物在 `08_figures/viz/`：

| 文件 | 内容 |
|---|---|
| `<作业>_complex_metrics.png/pdf/svg` | 四面板：复合物 RMSD、埋藏面积、配体埋藏比例（85% 参考线）、蛋白–配体质心距离 |
| `<作业>_snapshots_r*.png/pdf` | 0 / 中 / 末 三个时刻的结构快照，共用同一朝向 |
| `<作业>_bfactor_r*.png/pdf` | 逐残基 RMSF 映射到结构上的 B-factor 着色 + 色条 |
| `snapshots_r*.pml`、`bfactor_r*.pml` | **可手动编辑的 PyMOL 脚本**，改视角/配色后 `pymol xxx.pml` 重跑 |
| `*_snapshots_r*.pse`、`*_bfactor_r*.pse` | PyMOL 会话文件，File → Open 即可继续调 |

快照提取时会先把每帧按骨架叠合到首帧（`trjconv -fit rot+trans`），否则蛋白在盒子里自由转动，三个时刻的朝向没法比较。取景范围包含配体和它周围 12 Å 的蛋白，所以配体跑出口袋时也能看见它去了哪里。

埋藏 SASA 有两个口径：面积差（nm²）和配体被包住的百分比。后者文献里常用，>85% 通常认为埋得足够好。SASA 按 0.2 ns 抽帧算，100 ns 只有 500 帧，省时间且不影响趋势。

## 相互作用（interactions）

```bash
python scripts/mdpipe.py interactions --job Gal_1          # 全部重复
python scripts/mdpipe.py interactions --job Gal_1 --cutoff 0.4 --dt 0.5
python scripts/mdpipe.py interactions --job Gal_1 --no-ligplot
```

产物在 `08_figures/interactions/`：

| 文件 | 内容 |
|---|---|
| `contacts_r*.csv` | 逐残基接触占有率（全程 / 结合帧）、最小与平均距离、最近配体原子 |
| `contacts_r*.png/pdf` | 接触最频繁的 18 个残基（横向条形图） |
| `hbonds_r*.csv` | 逐对氢键：供体残基与原子、配体原子、占有率、平均距离与角度 |
| `hbonds_r*.png/pdf` | 占有率最高的 12 对氢键 |
| `clusters_r*.xvg/csv/png` | 配体构象聚类（GROMOS 法，cutoff 0.15 nm） |
| `interaction_map_r*.png/pdf/svg` | RDKit + matplotlib 的 2D 示意图，蓝虚线氢键、红点划线盐桥、灰点线疏水接触 |
| `ligplot_r*.png/pdf` | LigPlot+ 的经典 LIGPLOT 图（检测到 LigPlot+ 才出） |
| `ligplot_r*.pdb` | 接触残基最多的那一帧复合物，可直接拖进 LigPlot+ / PyMOL |
| `interactions_r*.json` | 上述结果的摘要 |

判据都是文献里的常规值：接触指任意重原子间距 ≤ 4.5 Å；氢键指 D···A ≤ 3.5 Å 且 D–H···A ≥ 120°；盐桥仅对羧基氧判定，葡萄糖的羟基不算。

占有率分两列：`occupancy` 是全程帧数占比，`occupancy_bound` 只在"配体与蛋白有接触"的帧里统计。弱结合体系（比如葡萄糖 20 ns 后离开口袋）看后者才有意义，图里画的也是后者。

**LigPlot+**：本机装在 `<LigPlot+目录>` 时脚本会自动调用它的 `hbplus.exe` + `ligplot.exe`（无需打开图形界面），再用 Ghostscript 转 PNG/PDF 并裁白边。换机器没装也不影响，会跳过这一步。注意命令行路径只能画氢键；LigPlot+ 图形界面里才有疏水接触的"睫毛"符号——需要的话打开 `ligplot_r*.pdb` 手动跑一次，或者直接看 `interaction_map_r*.png`（里面疏水接触按占有率画好了）。

## MM-PBSA 结合自由能

```bash
python scripts/mdpipe.py mmpbsa --job Gal_1 Glc_1            # 两个体系，全部重复
python scripts/mdpipe.py mmpbsa --job Gal_1 --frames 20 --last-ns 50
```

底层是 WSL 里的 gmx_MMPBSA（`/home/<用户名>/miniconda3/envs/mmpbsa`）。脚本替你处理了四件麻烦事：2019 版 gmx 读不了 2026 的 tpr、ParmEd 解析 `#include` 会读坏键参数、找不到 `amber14sb.ff`、以及原始轨迹里蛋白跨周期边界导致 BOND 项飙到 1e7（先 `trjconv -pbc mol -center`）。

取样窗口默认取轨迹**最后 20 ns**（`--last-ns`），这是文献里最常见的做法；`--last-ns 0` 表示全程平均——配体会解离的体系别用全程，否则等于在算"平均结合程度"。

产物：作业下 `08_mmpbsa_r*/mmpbsa.json`、`FINAL_RESULTS.dat`；多体系对比图 `pipeline/figures/MMPBSA_binding.png/pdf` 和同名 `.csv`。

> 逐残基分解（`--decomp`）目前不可用：这套 GROMACS→ParmEd 的拓扑会让 Amber 的内能项溢出成 `'*********'`，gmx_MMPBSA 直接报错。逐残基相互作用请用 `interactions` 命令的 `hbonds_r*.csv` / `contacts_r*.csv`。

## 自检与归档

每个作业分析完会自动跑一次自检，结论写进 `07_analysis*/selfcheck.json`，`status` 的「自检」列显示 `ok` / `warn` / `fail`。自检只用与 PBC 无关的量（Rg、SASA、接触原子数）交叉验证 RMSD，能抓出那类"RMSD 假涨、Rg 却正常"的轨迹问题。

体量最大的是能重算的副本。归档命令默认只列不删：

```bash
python scripts/mdpipe.py archive --all                    # 预演
python scripts/mdpipe.py archive --all --copy-to E:\分子动力学模拟   # 顺带校验 E 盘副本
python scripts/mdpipe.py archive --all --delete           # 真删
```

会清掉 `07_analysis*/md_center.xtc`、`md_prot.xtc`（`analysis` 可从 `md.xtc` 重算）和 `04_nvt*/nvt.xtc`、`05_npt*/npt.xtc`（平衡已结束）。缺 `md.xtc` 或 `md.tpr`、或者目标副本对不上的，一律跳过不删；每次删除都写 `pipeline/archive_<时间>.json` 清单。

## 完成通知与参数溯源

队列跑完会弹一条 Windows 桌面通知并响一声，失败时会说明是哪几个作业。关掉：

```json
{"notify": false, "notify_per_job": false}
```

`notify_per_job` 改 `true` 则每个作业完成/失败都提示一次（作业多时会吵）。单测通知：`python scripts/notify.py "标题" "内容"`。

每个作业在 `run` 开始时把运行环境写进 `job.json` 的 `provenance`，完成时补上全部 mdp 文件的哈希：

```bash
python scripts/mdpipe.py provenance --job Gal_1
```

包含 GROMACS 版本、AmberTools 版本、力场 `forcefield.itp` 的哈希、numpy/matplotlib 版本，以及 `mdpipe.py`、`selfcheck.py`、`make_figures.py`、`make_fel.py` 的代码指纹——以后改了分析口径，能一眼看出哪些结果是旧代码出的。老作业加 `--backfill` 补录，会标 `backfilled: true`，说明代码指纹是补录时的版本。

## 断点续跑

每一步成功后在作业目录里留一个 `.done_<步骤>` 标记。中断、关机、报错之后直接再跑 `run`，会从没做完的那一步接着走。失败的作业在 `job.json` 的 `error` 字段里留有最后一次命令的日志尾巴，完整输出在 `jobs/<作业名>/logs/`。

## GROMOS 54A7 + ATB

蛋白换成 GROMOS 时，配体参数基本只能从 ATB 网页下载。把导出的包直接指过去就行：

```bash
python scripts/mdpipe.py add --protein FABP5.pdb --ligand gal.pdbqt \
    --name Gal_G54 --forcefield gromos54a7_atb --ligand-atb ATB导出包/
```

**全原子导出必须配 ATB 的扩展力场 `gromos54a7_atb.ff`**（分子页右侧 Forcefield Files，免登录下载，解开丢进 `share/gromacs/top/`）：ATB 用的 `COpt`/`HOpt`/`OAlc`/`HS14` 在原生 54A7 里没有，缺了 grompp 直接报未定义类型。`atb_import.py --check` 会把缺的类型和该下哪个力场一起说出来。

`params` 步骤会调 `scripts/atb_import.py`：删掉与 `gromos54a7.ff` 重复的 `[ atomtypes ]`（不去重的话 grompp 直接报类型冲突）、统一 moleculetype 与残基名、核对 itp 与 gro 的原子数和顺序、缺 posre 就生成。mdp 也会自动换成 `mdp_templates/gromos/`（反应场 + 1.4 nm + `constraints = all-bonds`），不会套 amber 那套 PME + 1.0 nm。

GROMOS 作业的 grompp 会固定报一条"原本按双截断拟合"的警告，这条改 mdp 消除不掉，流水线只对 GROMOS 把 `-maxwarn` 放宽到 1。

## CHARMM-GUI

CHARMM-GUI 是网页服务，本流水线不代登录、不代提交。做法是在浏览器里走完向导，下载输入包，再在本机跑：

1. PDB Reader / Ligand Reader：清理结构、上传配体；
2. Solution Builder / Multicomponent Assembler：建溶剂、离子、盒子；
3. 下载输入包，解压后按包内 `README` 跑 `gmx grompp` / `gmx mdrun`。

注意 CHARMM-GUI 出的是 CHARMM36m + CGenFF，和本流水线的 amber14sb + GAFF2 不通用；默认的 GROMACS 安装通常也没有 CHARMM36m。详见手册 §16。

## 出图

每个作业跑完分析，会自动出图到 `jobs/<作业名>/08_figures/`：

```
<作业名>_MD_analysis.png / .pdf / .svg    六联图：骨架RMSD、配体RMSD、Rg、SASA、氢键、RMSF
<作业名>_MD统计.csv                        全程与后 50 ns 的均值±SD
```

一次重复红线、二次重复蓝线，叠在同一面板里。轨迹 ≥ 10 ns 时再多一步 PCA，出 `_PCA_FEL`（特征值谱、投影随时间、2D 景观）和 `_FEL_3D`（3D 曲面 + 地面映射）。想重出或补算老作业：

```bash
python scripts/mdpipe.py figures --job Gal_1      # 单体系
python scripts/mdpipe.py figures --all            # 队列里全部
python scripts/mdpipe.py figures --job Gal_1 Glc_1 Man_1   # 对比图 → pipeline/figures/
python scripts/mdpipe.py figures --job A B --labels "β-D-Gal" "D-GalA"   # 自定义图例
python scripts/mdpipe.py fel --job Gal_1 --force          # 短轨迹也强行算 FEL
```

界面里在作业列表选中作业（Ctrl 多选），点「出图」。分析口径（`md_prot.xtc` 即 `-pbc mol` 重整后的 Protein_LIG 轨迹 + 首帧参考、氢键取 `hbond` 的第 2 列）见手册 §20.7。

`figures` 只重画、不动 xvg。分析代码改过之后（或者作业是队列启动前登记的、跑的是内存里的旧代码），要重算分析本身：

```bash
python scripts/mdpipe.py analysis --job Gal_1 --replicates 1 2
```

## 测试时长

第一次用建议先跑一个短作业确认环境正常：

```bash
python scripts/mdpipe.py run --nvt-ps 50 --npt-ps 100 --ns 0.1
```

确认 EM 最大力 < 1000、NVT 温度 ≈ 300 K、NPT 密度 1020–1040 kg/m³ 之后，再按真实时长跑。

## 文件

| 路径 | 内容 |
|---|---|
| `config.json` | gmx 路径、WSL 发行版、AmberTools 路径、默认时长 |
| `queue.json` | 作业清单 |
| `jobs/<作业名>/` | 每个作业的完整 MD 目录（结构同手册） |
| `logs/`（在作业目录内） | 每步的完整终端输出 |

本目录下曾用 FABP5 + 半乳糖的两个位姿跑通全流程（结果存于 `review_checks/pipeline_smoke_summary.csv`），测试作业已清理。
