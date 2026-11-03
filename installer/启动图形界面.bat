@echo off
chcp 65001 >nul
rem 启动 mdpipe 图形界面。关掉这个黑窗口不会打断后台任务以外的任何东西。
setlocal
cd /d "%~dp0"

set "GUI=%~dp0scripts\mdpipe_gui.py"
set "PYW="
where pythonw >nul 2>nul && set "PYW=pythonw"
if not defined PYW (
    if exist "C:\Program Files\Python314\pythonw.exe" set "PYW=C:\Program Files\Python314\pythonw.exe"
)

if defined PYW (
    start "" "%PYW%" "%GUI%"
    exit /b 0
)

where python >nul 2>nul
if %errorlevel%==0 (
    python "%GUI%"
    pause
    exit /b 0
)

echo [错误] 找不到 Python。请先安装 Python 3.9 或更新版本。
pause
exit /b 1
