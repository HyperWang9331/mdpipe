#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""队列完成提示：Windows 桌面通知 + 蜂鸣。

只在 Windows 上弹通知（用 PowerShell 的 NotifyIcon，不依赖任何第三方包）；
别的平台上退化成蜂鸣和终端输出。所有失败都被吞掉——通知发不出去不能影响
MD 队列本身。

单独测试：
    python scripts/notify.py "标题" "内容"
"""

from __future__ import annotations

import base64
import shutil
import subprocess
import sys

PS_TEMPLATE = r"""
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing
$n = New-Object System.Windows.Forms.NotifyIcon
$n.Icon = [System.Drawing.SystemIcons]::{icon}
$n.Visible = $true
$n.ShowBalloonTip({ms}, '{title}', '{body}', '{kind}')
Start-Sleep -Milliseconds {hold}
$n.Dispose()
"""


def _ps_escape(s: str) -> str:
    return s.replace("'", "''").replace("\r", " ").replace("\n", " ")


def desktop(title: str, body: str, level: str = "info", hold_ms: int = 6000) -> bool:
    """弹一条桌面通知。返回是否成功发出。"""
    if sys.platform != "win32":
        return False
    exe = shutil.which("powershell") or shutil.which("pwsh")
    if not exe:
        return False
    icon, kind = {
        "info": ("Information", "Info"),
        "warn": ("Warning", "Warning"),
        "fail": ("Error", "Error"),
    }.get(level, ("Information", "Info"))
    script = PS_TEMPLATE.format(
        icon=icon, kind=kind, title=_ps_escape(title)[:120],
        body=_ps_escape(body)[:400], ms=10000, hold=hold_ms)
    try:
        subprocess.run(
            [exe, "-NoProfile", "-NonInteractive", "-WindowStyle", "Hidden",
             "-EncodedCommand", base64.b64encode(script.encode("utf-16-le")).decode()],
            timeout=20, capture_output=True, check=False)
        return True
    except Exception:
        return False


def beep(level: str = "info") -> None:
    if sys.platform == "win32":
        try:
            import winsound
            if level == "fail":
                winsound.Beep(620, 220)
                winsound.Beep(440, 320)
            elif level == "warn":
                winsound.Beep(700, 250)
            else:
                winsound.Beep(880, 160)
                winsound.Beep(1180, 220)
            return
        except Exception:
            pass
    try:
        sys.stdout.write("\a")
        sys.stdout.flush()
    except Exception:
        pass


def notify(title: str, body: str = "", level: str = "info",
           sound: bool = True, desktop_toast: bool = True) -> None:
    """发通知；任何环节失败都静默跳过。"""
    if sound:
        beep(level)
    if desktop_toast:
        desktop(title, body, level)
    if not desktop_toast or sys.platform != "win32":
        line = f"\n*** {title}" + (f"：{body}" if body else "") + " ***"
        print(line, flush=True)


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    argv = list(sys.argv[1:] if argv is None else argv)
    title = argv[0] if argv else "分子动力学流水线"
    body = argv[1] if len(argv) > 1 else "测试通知"
    level = argv[2] if len(argv) > 2 else "info"
    notify(title, body, level)
    print("已发送（Windows 上应能看到右下角通知）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
