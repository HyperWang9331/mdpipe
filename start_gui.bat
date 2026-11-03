@echo off
rem Launch the MD pipeline GUI. Paths are resolved from this file's folder,
rem so it works no matter which directory the shortcut starts in.
cd /d "%~dp0"
set "GUI=%~dp0scripts\mdpipe_gui.py"
set "PYW=C:\Program Files\Python314\pythonw.exe"

if exist "%PYW%" (
    start "" "%PYW%" "%GUI%"
    exit /b 0
)
where pythonw >nul 2>nul
if %errorlevel%==0 (
    start "" pythonw "%GUI%"
    exit /b 0
)
where python >nul 2>nul
if %errorlevel%==0 (
    python "%GUI%"
    pause
    exit /b 0
)
echo Python not found. Install Python or edit PYW in this file.
pause
