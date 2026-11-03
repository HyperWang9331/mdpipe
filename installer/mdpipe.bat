@echo off
chcp 65001 >nul
rem 命令行入口：mdpipe doctor / run / status ...
rem 路径以本文件所在目录为准，从任何工作目录调用都行。
setlocal
cd /d "%~dp0"

set "PY="
where python >nul 2>nul && set "PY=python"
if not defined PY (
    if exist "C:\Program Files\Python314\python.exe" set "PY=C:\Program Files\Python314\python.exe"
)
if not defined PY (
    echo [错误] 找不到 Python。请先安装 Python 3.9 或更新版本，安装时勾选 Add to PATH。
    pause
    exit /b 1
)

"%PY%" "%~dp0scripts\mdpipe.py" %*
