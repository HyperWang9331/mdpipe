# 打包 MDpipe：便携 zip + Inno Setup 安装程序
#
#   powershell -ExecutionPolicy Bypass -File installer\build.ps1
#   powershell -ExecutionPolicy Bypass -File installer\build.ps1 -ZipOnly
#
# 产物在 dist\ 下：
#   MDpipe-portable.zip       解压即用（目标机需自备 Python 3.9+）
#   MDpipe-Setup-1.7.exe      安装程序（开始菜单、卸载、加入 PATH 可选）

param(
    [switch]$ZipOnly,
    [string]$Iscc = ""
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
$dist = Join-Path $root "dist"
$stage = Join-Path $dist "MDpipe-portable"

Write-Host "== 1/3 收集文件 ==" -ForegroundColor Cyan
if (Test-Path $stage) { Remove-Item -Recurse -Force $stage }
New-Item -ItemType Directory -Force -Path $stage | Out-Null

Copy-Item (Join-Path $root "scripts")            (Join-Path $stage "scripts") -Recurse
Copy-Item (Join-Path $root "mdp_templates")      (Join-Path $stage "mdp_templates") -Recurse
if (Test-Path (Join-Path $root "tests")) {
    Copy-Item (Join-Path $root "tests") (Join-Path $stage "tests") -Recurse
}

New-Item -ItemType Directory -Force -Path (Join-Path $stage "pipeline") | Out-Null
foreach ($d in "inbox", "jobs", "figures", "splits") {
    New-Item -ItemType Directory -Force -Path (Join-Path $stage "pipeline\$d") | Out-Null
}
Copy-Item (Join-Path $root "pipeline\config.json") (Join-Path $stage "pipeline\config.json")
Copy-Item (Join-Path $root "pipeline\README.md")   (Join-Path $stage "pipeline\README.md")

New-Item -ItemType Directory -Force -Path (Join-Path $stage "skills\md-pipeline") | Out-Null
Copy-Item (Join-Path $root "skills\md-pipeline\SKILL.md") (Join-Path $stage "skills\md-pipeline")

New-Item -ItemType Directory -Force -Path (Join-Path $stage "docs") | Out-Null
Copy-Item (Join-Path $root "分子动力学操作手册.md") (Join-Path $stage "docs")
# 只带方法文档进安装包；README/LICENSE 留在仓库根目录
foreach ($doc in "安装与环境清单.md", "MD流程说明.md") {
    $src = Join-Path $root $doc
    if (Test-Path $src) { Copy-Item $src (Join-Path $stage "docs") }
}

Copy-Item (Join-Path $PSScriptRoot "mdpipe.bat")        (Join-Path $stage "mdpipe.bat")
Copy-Item (Join-Path $PSScriptRoot "启动图形界面.bat")  (Join-Path $stage "启动图形界面.bat")
Copy-Item (Join-Path $PSScriptRoot "环境自检.bat")      (Join-Path $stage "环境自检.bat")
Copy-Item (Join-Path $PSScriptRoot "安装说明.md")       (Join-Path $stage "docs\安装说明.md")
Copy-Item (Join-Path $PSScriptRoot "appicon.ico")       (Join-Path $stage "appicon.ico")

# 流程图与图标（文档里用相对路径 docs\pipeline\figures\ 引用）
New-Item -ItemType Directory -Force -Path (Join-Path $stage "docs\pipeline\figures") | Out-Null
Copy-Item (Join-Path $PSScriptRoot "appicon.png")       (Join-Path $stage "docs\pipeline\figures\appicon.png")
Copy-Item (Join-Path $root "pipeline\figures\installer_flow.png") (Join-Path $stage "docs\pipeline\figures")
Copy-Item (Join-Path $root "pipeline\figures\installer_flow.pdf") (Join-Path $stage "docs\pipeline\figures")

# .bat 必须是 CRLF：cmd 遇到 LF 会把中文 rem 行和下一行粘在一起，脚本直接报废
Get-ChildItem $stage -Recurse -Filter "*.bat" | ForEach-Object {
    $t = [System.IO.File]::ReadAllText($_.FullName)
    $t = ($t -replace "`r`n", "`n") -replace "`n", "`r`n"
    [System.IO.File]::WriteAllText($_.FullName, $t, (New-Object System.Text.UTF8Encoding($false)))
}

# 编译缓存和运行产物不进包
Get-ChildItem $stage -Recurse -Directory -Filter "__pycache__" |
    Remove-Item -Recurse -Force
Get-ChildItem $stage -Recurse -Include "*.pyc", "*.pyo" |
    Remove-Item -Force
Remove-Item (Join-Path $stage "pipeline\queue.json") -ErrorAction SilentlyContinue

Write-Host "== 2/3 打便携 zip ==" -ForegroundColor Cyan
$zip = Join-Path $dist "MDpipe-portable.zip"
if (Test-Path $zip) { Remove-Item -Force $zip }
Compress-Archive -Path (Join-Path $stage "*") -DestinationPath $zip -CompressionLevel Optimal
Write-Host ("   {0}  ({1:N1} MB)" -f $zip, ((Get-Item $zip).Length / 1MB))

if ($ZipOnly) {
    Write-Host "只打 zip，跳过安装程序（-ZipOnly）。" -ForegroundColor Yellow
    exit 0
}

Write-Host "== 3/3 编译安装程序 ==" -ForegroundColor Cyan
if (-not $Iscc) {
    $cand = @(
        (Join-Path $env:LOCALAPPDATA "Programs\Inno Setup 6\ISCC.exe"),
        "C:\Program Files (x86)\Inno Setup 6\ISCC.exe",
        "C:\Program Files\Inno Setup 6\ISCC.exe"
    )
    $Iscc = $cand | Where-Object { Test-Path $_ } | Select-Object -First 1
}
if (-not $Iscc -or -not (Test-Path $Iscc)) {
    Write-Host "找不到 ISCC.exe（Inno Setup 6）。" -ForegroundColor Red
    Write-Host "装好 Inno Setup 6 后重跑，或加 -Iscc <路径>；只想分发便携版就加 -ZipOnly。"
    exit 1
}

# .iss 里含中文，Inno 要求 UTF-8 带 BOM
$iss = Join-Path $PSScriptRoot "mdpipe.iss"
$text = [System.IO.File]::ReadAllText($iss)
if ($text -notmatch "\[Setup\]") {
    throw "mdpipe.iss 里没有 [Setup] 段，文件可能损坏"
}
[System.IO.File]::WriteAllText($iss, $text, (New-Object System.Text.UTF8Encoding($true)))

& $Iscc $iss
if ($LASTEXITCODE -ne 0) { throw "ISCC 编译失败（退出码 $LASTEXITCODE）" }

Get-ChildItem $dist -Filter "*.exe" | ForEach-Object {
    Write-Host ("   {0}  ({1:N1} MB)" -f $_.FullName, ($_.Length / 1MB))
}
Write-Host "完成。" -ForegroundColor Green
