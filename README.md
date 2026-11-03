# MDpipe：蛋白–配体分子动力学流水线

一套面向 Windows + WSL 的 GROMACS 蛋白–小分子 MD 自动化工具：从受体/配体文件出发，完成拓扑准备、能量最小化、NVT/NPT、生产、分析、MM-PBSA 与出图。命令行和图形界面功能一致。

> 仓库只包含代码、文档、测试和 mdp 模板，不包含任何轨迹、结构数据、ATB/CHARMM-GUI 下载内容或本机配置。

## 下载

Windows 用户不需要装 Python 就能先用起来，三种方式任选：

| 方式 | 适合 | 下载 |
|---|---|---|
| **安装包**（推荐） | 不想碰命令行 | [MDpipe-Setup-1.7.exe](https://github.com/HyperWang9331/mdpipe/releases/latest/download/MDpipe-Setup-1.7.exe)（3.0 MB） |
| **便携版** | 免安装、放 U 盘带走 | [MDpipe-portable.zip](https://github.com/HyperWang9331/mdpipe/releases/latest/download/MDpipe-portable.zip)（1.1 MB） |
| **源码** | 要读/改代码、自己打包 | `git clone`，见下方「快速开始」 |

全部历史版本在 [Releases](https://github.com/HyperWang9331/mdpipe/releases)。

**安装包做了什么**：装到 `%LOCALAPPDATA%\MDpipe`（**不需要管理员权限**），建开始菜单项和卸载项，可选加入 PATH、可选给 AI 助手装 `md-pipeline` 技能，装完自动跑一次环境自检。流程图见 [`pipeline/figures/installer_flow.png`](pipeline/figures/installer_flow.png)。

**安装包不做什么**：不装 Python、GROMACS、显卡驱动、WSL、conda，不下载任何运行时。装完先双击「环境自检」，缺什么它会逐项列出来。

> **Windows 提示“未知发布者”**：安装包没有代码签名（个人项目未购买证书），SmartScreen 会拦一下——点「更多信息」→「仍要运行」即可。

校验下载（PowerShell）：

```powershell
certutil -hashfile MDpipe-Setup-1.7.exe SHA256
# 应输出 5cb8d4c42c222c2c158124962934b9169d0f999578fd0f5989838088f8d49b0e
```

## 主要功能

- **准备**：受体体检与补链、配体位姿提取、GAFF2/AM1-BCC 参数化、拓扑组装、溶剂化与加离子
- **运行**：作业队列、断点续跑、双重复并行、GPU/CPU 切换、完成通知
- **分析**：RMSD、RMSF、Rg、SASA、氢键、接触、配体聚类、PCA/FEL
- **结合自由能**：WSL 内 gmx_MMPBSA（最后 N ns 取样），多体系对比
- **出图**：六联图、对比图、2D 相互作用图，PNG/PDF/SVG
- **工程化**：环境自检、路径自动探测、参数溯源、结果自检、归档、61 个单元测试
- **AI 协作**：`skills/md-pipeline/SKILL.md` 可供支持 Agent Skills 的工具读取

## 快速开始

前置：Windows 10/11、Python 3.9+、GROMACS（建议 CUDA 版）。MM-PBSA 与 ATB 参数化还需要 WSL2 + AmberTools。

```powershell
git clone https://github.com/HyperWang9331/mdpipe.git
cd mdpipe
python -m pip install numpy scipy matplotlib pandas

python scripts\mdpipe.py doctor                  # 环境自检，只读
python scripts\mdpipe.py config --detect --write # 探测并写入组件路径
python scripts\mdpipe.py gui                     # 图形界面
```

把蛋白 PDB 和配体位姿放进 `pipeline\inbox\<作业名>\`，或直接用：

```powershell
python scripts\mdpipe.py add --protein receptor.pdb --ligand pose.pdbqt --ns 100
python scripts\mdpipe.py run
python scripts\mdpipe.py status
```

完整说明见 [分子动力学操作手册.md](分子动力学操作手册.md)；组件安装见 [安装与环境清单.md](安装与环境清单.md)；流水线细节见 [pipeline/README.md](pipeline/README.md)。

## 自己打包（开发者）

`installer\` 下是 Inno Setup 源码，可生成按用户安装的安装程序（开始菜单、卸载、可选加入 PATH）。需要 [Inno Setup 6](https://jrsoftware.org/isdl.php)：

```powershell
powershell -ExecutionPolicy Bypass -File installer\build.ps1
```

安装时会自动运行 `mdpipe config --detect --write`，探测本机 GROMACS、WSL、AmberTools、gmx_MMPBSA 和 LigPlot+ 的路径。

## 第三方组件与数据

本仓库**不打包**下列软件或数据，需要用户自行安装/下载：

| 组件 | 用途 | 获取方式 |
|---|---|---|
| GROMACS（LGPL） | MD 引擎 | https://manual.gromacs.org |
| Python 包 numpy/scipy/matplotlib/pandas | 分析与出图 | PyPI |
| WSL2 + AmberTools | 配体参数化、MM-PBSA | conda-forge |
| gmx_MMPBSA | 结合自由能 | https://valdes-tresanco-ms.github.io/gmx_MMPBSA/ |
| PyMOL、LigPlot+、Ghostscript、OpenBabel | 可选出图/格式转换 | 各自官网，需接受其许可 |
| MODELLER | 可选补 loop | https://salilab.org/modeller/，需自行申请 license |
| ATB | GROMOS 配体参数 | https://atb.uq.edu.au，按网站条款下载与引用 |
| CHARMM-GUI | 膜蛋白/糖链/CHARMM36m 体系 | https://charmm-gui.org，在浏览器中操作 |

ATB 导出的拓扑、CHARMM-GUI 生成的输入包、PDB 结构文件都不在本仓库内；脚本只处理用户自己下载到本机的文件。

## 数据与结果

仓库不含任何模拟轨迹或结果。用本工具产生的数据请自行归档（例如 Zenodo、Figshare 或机构存储），并在论文中注明 DOI。`pipeline/jobs/` 已被 `.gitignore` 排除。

## 许可

代码采用 [MIT License](LICENSE)。文档中的方法描述可自由引用，请注明来源。

## 引用

如果本工具对你的工作有帮助，请引用仓库（`CITATION.cff`）以及实际使用的 GROMACS、AmberTools、gmx_MMPBSA、ATB 等组件。

## AI使用说明

本项目开发过程中使用了 ClawsGO 辅助代码实现与文档撰写。
