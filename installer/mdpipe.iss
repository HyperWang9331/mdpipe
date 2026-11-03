; MDpipe 安装脚本（Inno Setup 6）
; 用 build.ps1 编译：powershell -ExecutionPolicy Bypass -File installer\build.ps1
;
; 安装策略：
;   * 装到用户目录（{localappdata}），不需要管理员权限，安装后目录可写
;     —— 流水线要在 pipeline/jobs 下生成文件，Program Files 会因权限失败。
;   * 不打包 Python / GROMACS / WSL 环境，只装脚本和文档；缺什么由
;     “环境自检”逐项指出。GROMACS 是 LGPL、AmberTools 经 conda 分发，
;     这里不代发二进制。
;   * 卸载默认保留 pipeline\jobs 里的计算结果。

#define MyAppName "MDpipe"
#define MyAppVersion "1.7"
#define MyAppPublisher "MDpipe"
#define MyAppURL "https://github.com/HyperWang9331/mdpipe"
#define MyAppExeName "启动图形界面.bat"

[Setup]
AppId={{8C1B7A54-3D0E-4E2B-9C77-6A3B1D5F0E21}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
DefaultDirName={localappdata}\MDpipe
DefaultGroupName={#MyAppName}
SetupIconFile=appicon.ico
UninstallDisplayIcon={app}\appicon.ico
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir=..\dist
OutputBaseFilename=MDpipe-Setup-{#MyAppVersion}
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64compatible
ArchitecturesAllowed=x64compatible
UninstallDisplayName={#MyAppName}
; 关掉自动关闭的“安装完是否启动”复选框的默认行为，统一由 [Run] 控制
CloseApplications=no

[Languages]
Name: "zh"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "创建桌面快捷方式"; GroupDescription: "附加任务："; Flags: unchecked
Name: "addpath"; Description: "把 mdpipe 加入用户 PATH（可在任意目录运行 mdpipe 命令）"; GroupDescription: "附加任务："
Name: "a_skill"; Description: "为 AI 助手安装 md-pipeline 技能（写入 %USERPROFILE%\.clawsgo\skills，让 AI 会用 mdpipe）"; GroupDescription: "附加任务："

[Dirs]
Name: "{app}\pipeline\inbox"
Name: "{app}\pipeline\jobs"
Name: "{app}\pipeline\figures"
Name: "{app}\pipeline\splits"
Name: "{app}\logs"

[Files]
Source: "..\scripts\*.py"; DestDir: "{app}\scripts"; Excludes: "__pycache__\*"; Flags: ignoreversion recursesubdirs
Source: "..\scripts\*.sh"; DestDir: "{app}\scripts"; Flags: ignoreversion
Source: "..\pipeline\config.json"; DestDir: "{app}\pipeline"; Flags: onlyifdoesntexist uninsneveruninstall
Source: "..\pipeline\README.md"; DestDir: "{app}\pipeline"; Flags: ignoreversion
Source: "..\mdp_templates\*"; DestDir: "{app}\mdp_templates"; Flags: ignoreversion recursesubdirs
Source: "..\tests\*.py"; DestDir: "{app}\tests"; Flags: ignoreversion
Source: "..\分子动力学操作手册.md"; DestDir: "{app}\docs"; Flags: ignoreversion
Source: "..\安装与环境清单.md"; DestDir: "{app}\docs"; Flags: ignoreversion
Source: "..\安装包与流程优化规划.md"; DestDir: "{app}\docs"; Flags: ignoreversion skipifsourcedoesntexist
Source: "安装说明.md"; DestDir: "{app}\docs"; Flags: ignoreversion
; 流程图与文档里的相对链接保持一致：docs\pipeline\figures\installer_flow.png
Source: "..\pipeline\figures\installer_flow.png"; DestDir: "{app}\docs\pipeline\figures"; Flags: ignoreversion
Source: "..\pipeline\figures\installer_flow.pdf"; DestDir: "{app}\docs\pipeline\figures"; Flags: ignoreversion
Source: "..\skills\md-pipeline\SKILL.md"; DestDir: "{app}\skills\md-pipeline"; Flags: ignoreversion
Source: "..\skills\md-pipeline\SKILL.md"; DestDir: "{%USERPROFILE}\.clawsgo\skills\md-pipeline"; Flags: ignoreversion; Tasks: a_skill
Source: "appicon.ico"; DestDir: "{app}"; Flags: ignoreversion
Source: "appicon.png"; DestDir: "{app}\docs\pipeline\figures"; Flags: ignoreversion
Source: "mdpipe.bat"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#MyAppExeName}"; DestDir: "{app}"; Flags: ignoreversion
Source: "环境自检.bat"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\{#MyAppName}（图形界面）"; Filename: "{app}\{#MyAppExeName}"; WorkingDir: "{app}"; IconFilename: "{app}\appicon.ico"
Name: "{group}\环境自检"; Filename: "{app}\环境自检.bat"; WorkingDir: "{app}"; IconFilename: "{app}\appicon.ico"
Name: "{group}\操作手册"; Filename: "{app}\docs\分子动力学操作手册.md"; IconFilename: "{app}\appicon.ico"
Name: "{group}\安装与环境清单"; Filename: "{app}\docs\安装与环境清单.md"; IconFilename: "{app}\appicon.ico"
Name: "{group}\AI 操作说明（skill）"; Filename: "{app}\skills\md-pipeline\SKILL.md"; IconFilename: "{app}\appicon.ico"
Name: "{group}\命令行（mdpipe）"; Filename: "{cmd}"; Parameters: "/k ""{app}\mdpipe.bat"" --help"; WorkingDir: "{app}"; IconFilename: "{app}\appicon.ico"
Name: "{group}\卸载 {#MyAppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; WorkingDir: "{app}"; Tasks: desktopicon; IconFilename: "{app}\appicon.ico"

[Run]
Filename: "{app}\mdpipe.bat"; Parameters: "config --detect --write"; Description: "自动探测本机组件路径（GROMACS / WSL / AmberTools）"; Flags: postinstall runhidden
Filename: "{app}\环境自检.bat"; Description: "现在运行环境自检（推荐，只读，不装任何东西）"; Flags: postinstall shellexec skipifsilent
Filename: "{app}\{#MyAppExeName}"; Description: "启动 {#MyAppName}"; Flags: postinstall nowait skipifsilent unchecked

[UninstallDelete]
; 只删安装时创建的空目录和缓存；pipeline\jobs 是用户数据，不碰
Type: filesandordirs; Name: "{app}\scripts\__pycache__"
Type: filesandordirs; Name: "{app}\logs"
Type: dirifempty; Name: "{app}\pipeline\inbox"
Type: dirifempty; Name: "{app}\pipeline\splits"
Type: dirifempty; Name: "{app}\pipeline\figures"

[Code]
function PythonPresent(): Boolean;
var
  Names: TArrayOfString;
  I: Integer;
begin
  Result := False;
  if FileExists('C:\Program Files\Python314\python.exe') then
    Result := True;
  if (not Result) and FileExists('C:\Program Files\Python313\python.exe') then
    Result := True;
  if RegGetSubkeyNames(HKCU, 'SOFTWARE\Python\PythonCore', Names) then
    for I := 0 to GetArrayLength(Names) - 1 do
      if Pos('3.', Names[I]) = 1 then
        Result := True;
  if (not Result) and RegGetSubkeyNames(HKLM, 'SOFTWARE\Python\PythonCore', Names) then
    for I := 0 to GetArrayLength(Names) - 1 do
      if Pos('3.', Names[I]) = 1 then
        Result := True;
end;

procedure AddToUserPath();
var
  OldPath: String;
begin
  if not RegQueryStringValue(HKCU, 'Environment', 'Path', OldPath) then
    OldPath := '';
  if Pos(Lowercase(ExpandConstant('{app}')), Lowercase(OldPath)) = 0 then
  begin
    if (OldPath <> '') and (OldPath[Length(OldPath)] <> ';') then
      OldPath := OldPath + ';';
    RegWriteExpandStringValue(HKCU, 'Environment', 'Path', OldPath + ExpandConstant('{app}'));
  end;
end;

function InitializeSetup(): Boolean;
begin
  Result := True;
  if not PythonPresent() then
  begin
    if MsgBox('没有检测到 Python 3.9 或更新版本。' + #13#10#13#10 +
              '流水线本身是 Python 脚本，必须先装 Python。安装时请勾选 ' +
              '「Add python.exe to PATH」。装好后再运行本安装程序即可。' + #13#10#13#10 +
              '现在仍然继续安装（只复制文件，之后自己补 Python）吗？',
              mbConfirmation, MB_YESNO) = IDNO then
      Result := False;
  end;
end;

procedure CurStepChanged(CurStep: TSetupStep);
begin
  if (CurStep = ssPostInstall) and WizardIsTaskSelected('addpath') then
    AddToUserPath();
end;
