#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mdview.py —— 用 tkinter 在程序内渲染 Markdown，不跳外部软件。

只依赖标准库；装了 Pillow 时顺带把本地图片画出来，没装就显示占位文字。

    import mdview
    mdview.show(root, Path("docs/操作手册.md"))

支持的语法：标题、粗体/斜体、行内代码、围栏代码块、无序/有序列表、
引用、分隔线、表格、链接（http 用浏览器打开，本地 .md 在同一窗口打开）、
本地图片（PNG/JPG/GIF）。
"""

from __future__ import annotations

import os
import re
import unicodedata
import webbrowser
from pathlib import Path
from tkinter import Text, Tk, Toplevel, ttk
from tkinter import font as tkfont

# 行内记号：先代码，再图片、链接，最后粗体/斜体
TOKEN_RE = re.compile(
    r"(?P<code>`[^`\n]+`)"
    r"|(?P<image>!\[(?P<ialt>[^\]]*)\]\((?P<isrc>[^)\s]+)(?:\s+\"[^\"]*\")?\))"
    r"|(?P<link>\[(?P<ltext>[^\]]+)\]\((?P<lurl>[^)\s]+)(?:\s+\"[^\"]*\")?\))"
    r"|(?P<bold>\*\*(?P<btext>[^*\n]+)\*\*)"
    r"|(?P<italic>(?<!\*)\*(?P<itext>[^*\n]+)\*(?!\*))"
)

HEADING_RE = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
LIST_RE = re.compile(r"^(\s*)([-*+]|\d+[.)])\s+(.*)$")
HR_RE = re.compile(r"^\s*([-*_])\s*(?:\1\s*){2,}$")
FENCE_RE = re.compile(r"^\s*```")
TABLE_SEP_RE = re.compile(r"^\s*\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)+\|?\s*$")


def display_width(s: str) -> int:
    """按显示宽度算：中日韩全角字符占 2 列。"""
    return sum(2 if unicodedata.east_asian_width(ch) in "WF" else 1 for ch in s)


def pad(s: str, width: int) -> str:
    return s + " " * max(0, width - display_width(s))


def plain(s: str) -> str:
    """去掉行内记号，用于表格单元格对齐。"""
    def repl(m):
        if m.group("code"):
            return m.group("code")[1:-1]
        if m.group("image"):
            return m.group("ialt") or ""
        if m.group("link"):
            return m.group("ltext")
        if m.group("bold"):
            return m.group("btext")
        if m.group("italic"):
            return m.group("itext")
        return ""
    return TOKEN_RE.sub(repl, s)


# 折行时优先在这些字符之后断开，避免把 Ghostscript 拆成 G / hostscript
BREAK_AFTER = " 、，,；;：:/+"


def wrap_cell(s: str, width: int) -> list[str]:
    """按显示宽度折行；优先在空格/顿号/逗号处断开，实在没有才逐字断。"""
    width = max(1, width)
    out: list[str] = []
    cur = ""
    for ch in s:
        if display_width(cur) + display_width(ch) > width and cur:
            cut = max(cur.rfind(c) for c in BREAK_AFTER)
            if cut >= max(2, len(cur) // 3):
                out.append(cur[:cut + 1].rstrip())
                cur = cur[cut + 1:].lstrip()
            else:
                out.append(cur)
                cur = ""
        cur += ch
    out.append(cur)
    return out or [""]


def fit_widths(widths: list[int], total: int = 96, floor: int = 10) -> list[int]:
    """列太多/太宽时，从最宽的列开始削，保证整表能放进窗口。"""
    widths = list(widths)
    while sum(widths) + 2 * (len(widths) - 1) > total:
        i = widths.index(max(widths))
        if widths[i] <= floor:
            break
        widths[i] -= 1
    return widths


def split_row(line: str) -> list[str]:
    """把 | a | b | 拆成 ['a', 'b']，保留单元格里的转义竖线。"""
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|"):
        s = s[:-1]
    cells = re.split(r"(?<!\\)\|", s)
    return [c.strip().replace("\\|", "|") for c in cells]


def inline_tokens(s: str) -> list[tuple[str, str, str]]:
    """切成 (类型, 文本, 附加信息) 的序列；类型为 text/code/bold/italic/link/image。"""
    out: list[tuple[str, str, str]] = []
    pos = 0
    for m in TOKEN_RE.finditer(s):
        if m.start() > pos:
            out.append(("text", s[pos:m.start()], ""))
        if m.group("code"):
            out.append(("code", m.group("code")[1:-1], ""))
        elif m.group("image"):
            out.append(("image", m.group("ialt") or "", m.group("isrc")))
        elif m.group("link"):
            out.append(("link", m.group("ltext"), m.group("lurl")))
        elif m.group("bold"):
            out.append(("bold", m.group("btext"), ""))
        elif m.group("italic"):
            out.append(("italic", m.group("itext"), ""))
        pos = m.end()
    if pos < len(s):
        out.append(("text", s[pos:], ""))
    return out


def _pick_font(root) -> str:
    fams = set(tkfont.families(root))
    for f in ("Microsoft YaHei UI", "Microsoft YaHei", "Segoe UI", "SimHei"):
        if f in fams:
            return f
    return "TkDefaultFont"


class MarkdownWindow(Toplevel):
    """一个带滚动条的只读 Markdown 窗口。"""

    def __init__(self, master, path: Path):
        super().__init__(master)
        self.path = Path(path)
        self.title(self.path.name)
        self.geometry("980x760")
        self.minsize(620, 420)

        self.ui_font = _pick_font(self)
        self._images: list = []

        bar = ttk.Frame(self, padding=(8, 4))
        bar.pack(fill="x")
        ttk.Label(bar, text=str(self.path), foreground="#6B7785").pack(side="left")

        body = ttk.Frame(self)
        body.pack(fill="both", expand=True)
        self.text = Text(body, wrap="char", padx=28, pady=18,
                         relief="flat", borderwidth=0, cursor="arrow")
        sb = ttk.Scrollbar(body, command=self.text.yview)
        self.text.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        self.text.pack(fill="both", expand=True)

        self._config_tags()
        self.update_idletasks()          # 先拿到真实像素宽度，表格才知道能排多宽
        self._render(self.path.read_text(encoding="utf-8", errors="replace"))
        self.text.configure(state="disabled")

    # ---------------- 样式 ----------------
    def _config_tags(self):
        t = self.text
        base = tkfont.Font(family=self.ui_font, size=10)
        mono = "Consolas" if "Consolas" in tkfont.families(self) else "Courier New"
        mono_f = tkfont.Font(family=mono, size=9)

        # 等宽字体一个“显示单位”的像素宽度，用来把窗口宽度换算成列宽
        mono_9 = tkfont.Font(family=mono, size=9)
        self._unit_px = max(1, mono_9.measure("M"), mono_9.measure("中") / 2)
        t.configure(font=base, background="#FFFFFF", foreground="#24303D",
                    spacing1=2, spacing3=4)
        t.tag_configure("h1", font=(self.ui_font, 20, "bold"), spacing1=16,
                        spacing3=8, foreground="#17324F")
        t.tag_configure("h2", font=(self.ui_font, 15, "bold"), spacing1=14,
                        spacing3=6, foreground="#17324F")
        t.tag_configure("h3", font=(self.ui_font, 12, "bold"), spacing1=10,
                        spacing3=4)
        t.tag_configure("h4", font=(self.ui_font, 10, "bold"), spacing1=8)
        t.tag_configure("bold", font=(self.ui_font, 10, "bold"))
        t.tag_configure("italic", font=(self.ui_font, 10, "italic"))
        t.tag_configure("code", font=mono_f, background="#F1F3F6",
                        foreground="#8C3B2E")
        t.tag_configure("pre", font=mono_f, background="#F6F7F9",
                        lmargin1=18, lmargin2=18, rmargin=18,
                        spacing1=6, spacing3=6)
        t.tag_configure("quote", foreground="#5E6B7A", lmargin1=22,
                        lmargin2=22, background="#F7F9FB")
        t.tag_configure("li", lmargin1=24, lmargin2=42)
        t.tag_configure("hr", foreground="#C9D1DA", spacing1=8, spacing3=8)
        t.tag_configure("table", font=mono_f, lmargin1=12, lmargin2=12,
                        spacing1=0, spacing3=0)
        t.tag_configure("th", font=(mono, 9, "bold"), lmargin1=12, lmargin2=12,
                        spacing1=0, spacing3=0)
        t.tag_configure("muted", foreground="#6B7785")
        self._link_tag = 0

    # ---------------- 渲染 ----------------
    def _render(self, md: str):
        lines = md.splitlines()
        i, n = 0, len(lines)
        while i < n:
            line = lines[i]
            stripped = line.strip()

            if FENCE_RE.match(line):
                i = self._block_code(lines, i)
                continue

            if stripped.startswith("|") and i + 1 < n and TABLE_SEP_RE.match(lines[i + 1]):
                i = self._block_table(lines, i)
                continue

            m = HEADING_RE.match(stripped)
            if m:
                level = min(len(m.group(1)), 4)
                self._inline(m.group(2), f"h{level}")
                self.text.insert("end", "\n", f"h{level}")
                i += 1
                continue

            if HR_RE.match(line):
                self.text.insert("end", "─" * 58 + "\n", "hr")
                i += 1
                continue

            if stripped.startswith(">"):
                self._block_quote(lines, i)
                i = self._quote_end(lines, i)
                continue

            m = LIST_RE.match(line)
            if m:
                self._block_list(lines, i)
                i = self._list_end(lines, i)
                continue

            if not stripped:
                i += 1
                continue

            # 普通段落：吃到空行或下一个块级记号
            para = []
            while i < n:
                s = lines[i].strip()
                if (not s or FENCE_RE.match(lines[i]) or HEADING_RE.match(s)
                        or HR_RE.match(lines[i]) or s.startswith(">")
                        or LIST_RE.match(lines[i])):
                    break
                para.append(s)
                i += 1
            self._inline(" ".join(para), "")
            self.text.insert("end", "\n")
            self.text.insert("end", "\n")

    def _inline(self, s: str, tag: str):
        for kind, text, extra in inline_tokens(s):
            if kind == "text":
                self.text.insert("end", text, tag or "body")
            elif kind in ("bold", "italic"):
                # 标题已经够粗够大，再加 bold/italic 会把标题字号覆盖掉
                use = tag if tag.startswith("h") else ((tag, kind) if tag else kind)
                self.text.insert("end", text, use)
            elif kind == "code":
                self.text.insert("end", text, ("code", tag))
            elif kind == "link":
                self._insert_link(text, extra)
            elif kind == "image":
                self._insert_image(text, extra)

    def _insert_link(self, label: str, url: str):
        self._link_tag += 1
        name = f"link{self._link_tag}"
        self.text.tag_configure(name, foreground="#2A6EA8", underline=True)
        self.text.insert("end", label, name)
        self.text.tag_bind(name, "<Button-1>", lambda e, u=url: self._open_link(u))
        self.text.tag_bind(name, "<Enter>",
                           lambda e: self.text.configure(cursor="hand2"))
        self.text.tag_bind(name, "<Leave>",
                           lambda e: self.text.configure(cursor="arrow"))

    def _open_link(self, url: str):
        if url.startswith(("http://", "https://", "mailto:")):
            webbrowser.open(url)
            return
        target = (self.path.parent / url.split("#")[0]).resolve()
        if target.suffix.lower() == ".md" and target.is_file():
            MarkdownWindow(self.master, target)
        elif target.is_file():
            try:
                import os
                os.startfile(str(target))
            except Exception:
                pass

    def _insert_image(self, alt: str, src: str):
        if src.startswith(("http://", "https://")):
            self.text.insert("end", f"[图片] {alt or src}", "muted")
            return
        p = (self.path.parent / src).resolve()
        photo = None
        if p.is_file():
            try:
                from PIL import Image, ImageTk
                im = Image.open(p)
                w, h = im.size
                max_w = max(320, self._content_px() - 34)
                if w > max_w:
                    im = im.resize((max_w, max(1, int(h * max_w / w))),
                                   Image.LANCZOS)
                photo = ImageTk.PhotoImage(im)
            except Exception:
                photo = None
        if photo is None:
            self.text.insert("end", f"[图片] {alt or p.name}", "muted")
            return
        self._images.append(photo)
        self.text.insert("end", "\n")
        self.text.image_create("end", image=photo)
        self.text.insert("end", "\n")
        if alt:
            self.text.insert("end", alt + "\n", "muted")

    def _content_px(self) -> int:
        """正文可用像素宽度。窗口还没真正显示时退回 geometry() 里的宽度。"""
        w = self.text.winfo_width()
        if w <= 1:
            try:
                w = int(self.winfo_geometry().split("x")[0])
            except Exception:
                w = 980
        return max(320, w - 56)

    def _table_total(self) -> int:
        """当前窗口一行大约放得下多少个等宽显示单位。"""
        return max(40, int((self._content_px() - 30) / self._unit_px))

    def _block_code(self, lines: list[str], i: int) -> int:
        i += 1
        body = []
        while i < len(lines) and not FENCE_RE.match(lines[i]):
            body.append(lines[i])
            i += 1
        self.text.insert("end", "\n" + "\n".join(body) + "\n", "pre")
        return i + 1

    def _block_quote(self, lines: list[str], i: int):
        while i < len(lines) and lines[i].strip().startswith(">"):
            content = lines[i].strip()[1:].strip()
            self.text.insert("end", "│ ", "quote")
            self._inline(content, "quote")
            self.text.insert("end", "\n", "quote")
            i += 1

    def _quote_end(self, lines: list[str], i: int) -> int:
        while i < len(lines) and lines[i].strip().startswith(">"):
            i += 1
        return i

    def _block_list(self, lines: list[str], i: int):
        while i < len(lines):
            m = LIST_RE.match(lines[i])
            if not m:
                break
            indent, marker, content = m.groups()
            bullet = "•" if not marker[0].isdigit() else marker.rstrip(".)")
            prefix = "    " * (len(indent) // 2)
            self.text.insert("end", f"{prefix}{bullet}  ", "li")
            self._inline(content, "li")
            self.text.insert("end", "\n", "li")
            i += 1

    def _list_end(self, lines: list[str], i: int) -> int:
        while i < len(lines) and LIST_RE.match(lines[i]):
            i += 1
        return i

    def _block_table(self, lines: list[str], i: int) -> int:
        header = [plain(c) for c in split_row(lines[i])]
        i += 2
        rows = []
        while i < len(lines) and lines[i].strip().startswith("|"):
            rows.append([plain(c) for c in split_row(lines[i])])
            i += 1
        ncol = max([len(header)] + [len(r) for r in rows])
        header += [""] * (ncol - len(header))
        rows = [r + [""] * (ncol - len(r)) for r in rows]

        widths = [max(display_width(header[c]),
                      *[display_width(r[c]) for r in rows]) for c in range(ncol)]
        widths = fit_widths(widths, total=self._table_total(), floor=8)

        def emit(cells, tag):
            wrapped = [wrap_cell(cells[c], widths[c]) for c in range(ncol)]
            for k in range(max(len(w) for w in wrapped)):
                line = "  ".join(
                    pad(w[k] if k < len(w) else "", widths[c])
                    for c, w in enumerate(wrapped))
                self.text.insert("end", line.rstrip() + "\n", tag)

        self.text.insert("end", "\n")
        emit(header, "th")
        self.text.insert("end", "  ".join("─" * w for w in widths) + "\n", "table")
        for r in rows:
            emit(r, "table")
        self.text.insert("end", "\n")
        return i


def show(master, path: Path) -> MarkdownWindow:
    return MarkdownWindow(master, path)


def main(argv: list[str] | None = None) -> int:
    import sys
    args = sys.argv[1:] if argv is None else argv
    if not args:
        print("用法: python scripts/mdview.py <文件.md>")
        return 2
    root = Tk()
    root.withdraw()
    MarkdownWindow(root, Path(args[0]))
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
