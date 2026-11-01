"""
Helpers for building the manuscript and supplement as JSON for build_docx.js.

Inline markup in text strings:
    **bold**   *italic*   ^{superscript}   _{subscript}   `monospace`
Citations: [@key] or [@key1,@key2] -> [n] / [n,m] / [n-m], numbered by order
of first appearance (the reference list follows that order).
"""
import json
import re

TOKEN = re.compile(r"(\*\*.+?\*\*|\*[^*\s][^*]*?\*|\^\{.+?\}|_\{.+?\}|`.+?`)")


def parse_inline(s):
    out = []
    for part in TOKEN.split(s):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            for r in parse_inline(part[2:-2]):
                r["b"] = True
                out.append(r)
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            for r in parse_inline(part[1:-1]):
                r["i"] = True
                out.append(r)
        elif part.startswith("^{"):
            out.append({"t": part[2:-1], "sup": True})
        elif part.startswith("_{"):
            out.append({"t": part[2:-1], "sub": True})
        elif part.startswith("`"):
            out.append({"t": part[1:-1], "mono": True})
        else:
            out.append({"t": part})
    return out


class Citer:
    def __init__(self, refs):
        self.refs = refs            # key -> formatted reference string
        self.order = []

    def num(self, key):
        if key not in self.refs:
            raise KeyError(f"unknown reference key: {key}")
        if key not in self.order:
            self.order.append(key)
        return self.order.index(key) + 1

    def sub(self, text):
        def repl(m):
            keys = [k.strip().lstrip("@") for k in m.group(1).split(",")]
            nums = sorted(set(self.num(k) for k in keys))
            # compress runs of 3+ consecutive numbers
            parts, i = [], 0
            while i < len(nums):
                j = i
                while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                    j += 1
                if j - i >= 2:
                    parts.append(f"{nums[i]}–{nums[j]}")
                else:
                    parts.extend(str(n) for n in nums[i:j + 1])
                i = j + 1
            return "[" + ",".join(parts) + "]"
        return re.sub(r"\[(@[^\]]+)\]", repl, text)

    def reference_list(self):
        return [f"{i + 1}. {self.refs[k]}" for i, k in enumerate(self.order)]


class Doc:
    def __init__(self, citer=None, **meta):
        self.blocks = []
        self.meta = meta
        self.c = citer
        self.md = []

    def _t(self, text):
        text = self.c.sub(text) if self.c else text
        text = re.sub(r"(?<![\w\)\]/])-(?=\d)", "\u2212", text)   # unicode minus for numbers
        text = text.replace("= < ", "< ").replace("p = <", "p <")
        return text

    def h(self, level, text):
        self.blocks.append({"type": f"h{level}", "text": text})
        self.md.append("#" * level + " " + text)

    def p(self, text, style="body"):
        t = self._t(text)
        self.blocks.append({"type": "para", "style": style, "runs": parse_inline(t)})
        self.md.append(t)

    def table(self, caption, header, rows, widths=None, notes=None, font_size=8.5, align=None,
              shade_rows=None):
        cap = self._t(caption)
        self.blocks.append({"type": "para", "style": "tablecaption", "runs": parse_inline(cap)})
        self.md.append(cap)
        hdr = [[parse_inline(self._t(str(c))) for c in r] for r in header]
        body = [[parse_inline(self._t(str(c))) for c in r] for r in rows]
        self.blocks.append({"type": "table", "header": hdr, "rows": body, "widths": widths,
                            "font_size": font_size, "align": align, "shade_rows": shade_rows or []})
        for r in header + rows:
            self.md.append(" | ".join(str(c) for c in r))
        if notes:
            n = self._t(notes)
            self.blocks.append({"type": "para", "style": "tablenote", "runs": parse_inline(n)})
            self.md.append(n)

    def figure(self, path, caption, width_in=6.7, alt=""):
        self.blocks.append({"type": "figure", "path": path, "width_in": width_in, "alt": alt})
        cap = self._t(caption)
        self.blocks.append({"type": "para", "style": "caption", "runs": parse_inline(cap)})
        self.md.append(cap)

    def pagebreak(self):
        self.blocks.append({"type": "pagebreak"})

    def section(self, landscape=False):
        self.blocks.append({"type": "section", "landscape": landscape})

    def refs(self):
        for r in self.c.reference_list():
            self.blocks.append({"type": "para", "style": "ref", "runs": parse_inline(r)})
            self.md.append(r)

    def save(self, json_path, md_path):
        json.dump({"meta": self.meta, "blocks": self.blocks}, open(json_path, "w"), ensure_ascii=False)
        open(md_path, "w").write("\n\n".join(self.md))
