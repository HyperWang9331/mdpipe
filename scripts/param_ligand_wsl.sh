#!/bin/bash
# ============================================================================
# param_ligand_wsl.sh —— 在 WSL 内用 GAFF2 + AM1-BCC 生成配体 GROMACS 拓扑
#
# 用法 (在 WSL 里执行):
#   bash /path/to/param_ligand_wsl.sh 输入.mol2 残基名 净电荷
# 例:
#   bash param_ligand_wsl.sh pose_h.mol2 LIG 0
#
# 产物: <残基名>.itp, posre_<残基名>.itp, <残基名>.gro
# ============================================================================
set -euo pipefail

# AmberTools 位置：优先用环境变量 AMBERHOME，其次找常见的 conda 环境
if [ -z "${AMBERHOME:-}" ]; then
    for cand in "$HOME/miniconda3/envs/amber" "$HOME/miniforge3/envs/amber" \
                "$HOME/anaconda3/envs/amber" /opt/amber /usr/local/amber; do
        if [ -x "$cand/bin/antechamber" ]; then
            AMBERHOME="$cand"
            break
        fi
    done
fi
if [ -z "${AMBERHOME:-}" ] || [ ! -x "$AMBERHOME/bin/antechamber" ]; then
    echo "找不到 AmberTools。请先设置 AMBERHOME，例如：" >&2
    echo "  AMBERHOME=$HOME/miniconda3/envs/amber bash param_ligand_wsl.sh pose.mol2 LIG 0" >&2
    exit 1
fi
export AMBERHOME
export PATH="$AMBERHOME/bin:$PATH"

MOL2=${1:?用法: param_ligand_wsl.sh ligand.mol2 RESNAME CHARGE}
RES=${2:-LIG}
CHG=${3:-0}

echo "=== 环境 ==="
which antechamber
which acpype

echo "=== 1/2 antechamber: AM1-BCC 电荷 + GAFF2 原子类型 ==="
antechamber -i "$MOL2" -fi mol2 -o "${RES}_bcc.mol2" -fo mol2 \
    -c bcc -nc "$CHG" -rn "$RES" -at gaff2 -pf y

echo "=== 2/2 acpype: 转 GROMACS 拓扑 ==="
acpype -i "${RES}_bcc.mol2" -b "$RES" -c user -a gaff2 -n "$CHG"

echo
echo "=== 产物 ==="
ls -1 "${RES}".itp "posre_${RES}.itp" "${RES}".gro 2>/dev/null || ls -1
