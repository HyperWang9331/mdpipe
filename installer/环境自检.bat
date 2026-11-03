@echo off
chcp 65001 >nul
rem 只读检测：GROMACS / GPU / WSL / AmberTools / 出图工具 / 凭据 / 磁盘。
rem 同时把报告写到 doctor_report.txt，方便发给别人看。
setlocal
cd /d "%~dp0"

set "PY="
where python >nul 2>nul && set "PY=python"
if not defined PY (
    if exist "C:\Program Files\Python314\python.exe" set "PY=C:\Program Files\Python314\python.exe"
)
if not defined PY (
    echo [错误] 找不到 Python。请先安装 Python 3.9 或更新版本。
    pause
    exit /b 1
)

"%PY%" "%~dp0scripts\doctor.py" --report "%~dp0doctor_report.txt" %*
pause
