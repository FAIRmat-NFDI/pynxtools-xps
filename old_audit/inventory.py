#!/usr/bin/env python3
"""
inventory.py - build a review ledger for one pynxtools plugin.

Usage (from the plugin's root, next to pyproject.toml):
    python inventory.py            # writes inventory.csv
    python inventory.py /path/to/plugin
Requires Python >= 3.11 (tomllib).
"""
import ast
import csv
import re
import sys
import tomllib
from pathlib import Path

PLUGIN = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
SKIP_DIRS = {"__pycache__", ".git", ".venv", "venv", "build", "dist", "node_modules"}

PATTERNS = {
    "unit":    r"units?\b|@units|ureg|pint|\beV\b|degree|kelvin|\bum\b|counts",
    "math":    r"\bnp\.\w+|numpy|scipy|\.mean\(|\.sum\(|\.std\(|round\(|flip|interp|cumsum|\* ?1e|/ ?1000",
    "default": r"\.get\([^)]*,[^)]*\)|fallback|default|\belse\b ['\"\d]",
    "mapping": r"_MAP\b|key_map|value_map|unit_map|_context",
    "todo":    r"TODO|FIXME|XXX|TBD|^\s*pass\s*$",
    "raises":  r"\braise\b|_logger\.(warning|error)|warnings\.warn",
}

COLUMNS = ["plugin", "file", "kind", "qualname", "line_start", "line_end", "loc",
           "auto_flags", "status", "finding_ids", "reviewer_notes"]
rows = []


def flags(src: str) -> str:
    return ";".join(k for k, p in PATTERNS.items() if re.search(p, src, re.M))


def skipped(path: Path) -> bool:
    return any(part in SKIP_DIRS for part in path.relative_to(PLUGIN).parts)


def add(file, kind, qualname, start="", end="", src=""):
    loc = (end - start + 1) if isinstance(start, int) and isinstance(end, int) else ""
    rows.append([PLUGIN.name, file, kind, qualname, start, end, loc,
                 flags(src) if src else "", "unreviewed", "", ""])


def assigned_names(node):
    """Names bound by Assign / AnnAssign (plain names only)."""
    if isinstance(node, ast.Assign):
        targets = node.targets
    elif isinstance(node, ast.AnnAssign):
        targets = [node.target]
    else:
        return []
    return [t.id for t in targets if isinstance(t, ast.Name)]


def walk(node, rel, lines, prefix="", in_class=False):
    for ch in ast.iter_child_nodes(node):
        if isinstance(ch, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            kind = "class" if isinstance(ch, ast.ClassDef) else (
                "method" if in_class else "function")
            q = f"{prefix}{ch.name}"
            seg = "\n".join(lines[ch.lineno - 1: ch.end_lineno])
            add(rel, kind, q, ch.lineno, ch.end_lineno, seg)
            if isinstance(ch, ast.ClassDef):
                walk(ch, rel, lines, q + ".", in_class=True)
            else:
                walk(ch, rel, lines, q + ".<locals>.", in_class=False)  # nested defs
        elif isinstance(ch, (ast.Assign, ast.AnnAssign)):
            # module-level constants/instances and class attributes
            if prefix == "" or (in_class and not prefix.endswith("<locals>.")):
                for name in assigned_names(ch):
                    seg = "\n".join(lines[ch.lineno - 1: ch.end_lineno])
                    add(rel, "class_attribute" if in_class else "module_constant",
                        f"{prefix}{name}", ch.lineno, ch.end_lineno, seg)


# 1. entry points and metadata from pyproject.toml
pp = PLUGIN / "pyproject.toml"
if pp.exists():
    data = tomllib.loads(pp.read_text())
    for group, items in data.get("project", {}).get("entry-points", {}).items():
        for name, target in items.items():
            add("pyproject.toml", "entry_point", f"{group}:{name} = {target}")
    for dep in data.get("project", {}).get("dependencies", []):
        add("pyproject.toml", "dependency", dep)
else:
    print("WARNING: no pyproject.toml found in", PLUGIN)

# 2. Python code under src/ (the package)
for f in sorted((PLUGIN / "src").rglob("*.py")):
    if skipped(f):
        continue
    rel = str(f.relative_to(PLUGIN))
    text = f.read_text(errors="ignore")
    try:
        tree = ast.parse(text)
    except SyntaxError as e:
        add(rel, "PARSE_ERROR", str(e))
        continue
    add(rel, "module", f.name, 1, len(text.splitlines()), text)  # module-level view
    walk(tree, rel, text.splitlines())

# 3. non-Python files that affect behaviour or define expectations
GROUPS = [
    ("src/**/*.json", "config"), ("src/**/*.yaml", "config"), ("src/**/*.yml", "config"),
    ("src/**/*.toml", "config"), ("tests/**/*.py", "test_file"),
    ("docs/**/*.md", "doc_page"), ("examples/**/*", "example_file"),
    ("tests/data/**/*", "test_data_file"),
]
seen = set()
for pattern, kind in GROUPS:
    for f in sorted(PLUGIN.glob(pattern)):
        if not f.is_file() or skipped(f) or f in seen or f.suffix == ".pyc":
            continue
        seen.add(f)
        rel = str(f.relative_to(PLUGIN))
        try:
            n = sum(1 for _ in open(f, "rb"))
        except OSError:
            n = ""
        add(rel, kind, f.name, 1 if n != "" else "", n, "")

with open(PLUGIN / "inventory.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(COLUMNS)
    w.writerows(rows)

# summary for completeness checks
from collections import Counter
counts = Counter(r[2] for r in rows)
print(f"{len(rows)} units written to {PLUGIN / 'inventory.csv'}")
for kind, n in sorted(counts.items()):
    print(f"  {kind:16} {n}")
