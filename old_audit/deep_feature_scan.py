#!/usr/bin/env python3
"""
unorthodox_scan.py - static scanner that flags suspicious behaviour in a pynxtools
plugin (or any Python package): hard-coded values, unusual units, silent
fallbacks, data-changing operations, swallowed exceptions, copy-pasted logic.

It does NOT prove a bug. Every finding is a pointer for manual review.

Usage (from the plugin root, next to pyproject.toml; Python >= 3.9):
    python unorthodox_scan.py
    python unorthodox_scan.py --inventory inventory.csv      # also annotate the ledger
    python unorthodox_scan.py --include-tests --min-severity low
    python unorthodox_scan.py --disable D01,D05,A03          # silence noisy rules

Outputs (written next to the script's working directory):
    unorthodox_findings.csv           one row per finding
    inventory_with_findings.csv  (only with --inventory) ledger + counts per unit
"""
from __future__ import annotations

import argparse
import ast
import csv
import difflib
import hashlib
import io
import json
import re
import sys
import tokenize
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

SEV = {"info": 0, "low": 1, "medium": 2, "high": 3}
SKIP_DIRS = {".git", "__pycache__", ".venv", "venv", "build", "dist",
             "node_modules", ".mypy_cache", ".ruff_cache", ".pytest_cache"}

RULES = {
    # units
    "U01": "unit string not recognised",
    "U02": "ambiguous unit string",
    "U03": "non-canonical unit alias",
    "U04": "non-SI unit (check NeXus type / pint)",
    "U05": "unit string with whitespace or empty",
    "U06": "unit label hard-coded in code (not derived from data)",
    "U07": "unit hard-coded in config",
    # numbers
    "N01": "scale factor / physical constant in arithmetic",
    "N02": "magic number in arithmetic",
    "N03": "numeric default parameter",
    "N04": "rounding to fixed decimals",
    "N05": "truncating cast / split",
    "N06": "silent truncation to shorter length",
    # data-changing operations
    "A01": "axis flip / reverse",
    "A02": "numerical transform (interpolation, filter, cumulative ...)",
    "A03": "reduction (sum / mean / std ...)",
    # defaults & fallbacks
    "D01": ".get() with literal default",
    "D02": "option default read from kwargs",
    "D03": "'or <literal>' fallback",
    "D04": "constant array created (possible fabricated data)",
    "D05": "hard-coded string comparison (case-sensitive)",
    # errors
    "E01": "bare except",
    "E02": "exception swallowed",
    "E03": "assert used for validation",
    "E04": "unimplemented stub",
    "E05": "TODO / FIXME / TBD marker",
    # time
    "T01": "timezone assumption",
    "T02": "sentinel date",
    "T03": "strptime / naive datetime",
    # consistency
    "K01": "duplicate dict key (last one silently wins)",
    "K02": "keys differing only by case/underscore/space",
    "K03": "diverged copies of a function",
    "K04": "identical copies of a function",
    # config / paths
    "C01": "NeXus path literal in Python code",
    "C02": "hard-coded string value in config",
    "C03": "hard-coded number in config",
    "C04": "duplicate key in JSON object",
    "C05": "config file does not parse",
    # misc
    "M01": "print() call",
}

# ----------------------------------------------------------------- units ----
SI_PREFIXES = {"", "y", "z", "a", "f", "p", "n", "u", "\u00b5", "\u03bc", "m",
               "c", "d", "k", "M", "G", "T", "P"}
PREFIXABLE = {"eV", "V", "A", "W", "m", "s", "Hz", "Pa", "J", "K", "N", "T", "F", "S"}
CANON_ATOMS = {"eV", "J", "K", "V", "A", "W", "m", "s", "Hz", "Pa", "min", "h",
               "degree", "radian", "rad", "counts", "count", "counts_per_second",
               "angstrom", "dimensionless", "percent", "%", "degC", "N", "T", "S",
               "F", "ohm", "1"}
AMBIGUOUS = {"C": "coulomb or Celsius?", "nU": "nano-units; voltage or current?",
             "d": "VAMAS 'd' (unspecified)", "u": "atomic mass unit or micro?",
             "a.u.": "arbitrary units", "au": "arbitrary / atomic units",
             "arb": "arbitrary units", "arb.": "arbitrary units",
             "norm": "normalised (dimensionless?)", "unitless": "dimensionless?",
             "units": "placeholder", "unit": "placeholder"}
ALIASES = {"deg": "degree", "degrees": "degree", "\u00b0": "degree",
           "cps": "counts_per_second", "c/s": "counts_per_second",
           "counts/s": "counts_per_second", "counts per second": "counts_per_second",
           "\u00b5m": "um", "\u03bcm": "um", "\u00c5": "angstrom", "A\u030a": "angstrom",
           "sec": "s", "secs": "s", "seconds": "s", "hr": "h", "hrs": "h",
           "volt": "V", "volts": "V", "ev": "eV", "EV": "eV", "Ev": "eV",
           "kelvin": "K", "kev": "keV", "KeV": "keV", "mev": "meV", "S-1": "1/s",
           "s-1": "1/s", "s^-1": "1/s", "Hertz": "Hz", "hz": "Hz", "ua": "uA"}
NON_SI = {"mbar", "bar", "torr", "Torr", "atm", "psi", "mmHg", "inch", "in", "mil",
          "fahrenheit", "degF"}


def _atom_ok(atom: str) -> bool:
    atom = re.sub(r"\^?[-+]?\d+$", "", atom.strip())
    if atom == "" or re.fullmatch(r"[-+]?\d+(\.\d+)?", atom):
        return True
    if atom in CANON_ATOMS:
        return True
    for base in PREFIXABLE:
        if atom.endswith(base) and atom[: -len(base)] in SI_PREFIXES:
            return True
    return False


def classify_unit(u: str):
    """Return list of (rule, severity, message) for a unit string."""
    out = []
    if u == "":
        return [("U05", "low", "empty unit string")]
    if u != u.strip():
        out.append(("U05", "medium", f"unit has surrounding whitespace: {u!r}"))
        u = u.strip()
    if u in AMBIGUOUS:
        out.append(("U02", "medium", f"ambiguous unit {u!r}: {AMBIGUOUS[u]}"))
    elif u in ALIASES:
        out.append(("U03", "low", f"non-canonical alias {u!r}; prefer {ALIASES[u]!r}"))
    elif u in NON_SI:
        out.append(("U04", "low", f"non-SI unit {u!r}"))
    else:
        atoms = [a for a in re.split(r"[/*\u00b7]", u) if a != ""]
        if not atoms or not all(_atom_ok(a) for a in atoms):
            out.append(("U01", "medium", f"unrecognised unit string {u!r}"))
    return out


# ------------------------------------------------------------- numbers -----
SCALE_VALUES = [1e3, 1e6, 1e9, 1e12, 1e-3, 1e-6, 1e-9, 1e-12, 60.0, 3600.0, 86400.0, 100.0]
PHYS_VALUES = {  # approximate, 1 % tolerance
    1.602176634e-19: "elementary charge / eV->J",
    6.62607015e-34: "Planck constant (J s)",
    299792458.0: "speed of light",
    12398.419843: "hc in eV*Angstrom",
    1239.84198: "hc in eV*nm",
    8.617333262e-5: "Boltzmann (eV/K)",
    273.15: "Celsius -> Kelvin",
    1.0e-10: "Angstrom -> m",
}
IGNORE_NUMBERS = {0, 1, -1, 2, 0.5}


def number_hint(v):
    av = abs(float(v))
    for s in SCALE_VALUES:
        if av == s:
            return "N01", "medium", f"possible unit/scale factor {v!r}"
    for p, label in PHYS_VALUES.items():
        if p != 0 and abs(av - p) / p < 0.01:
            return "N01", "medium", f"looks like a physical constant ({label}): {v!r}"
    return "N02", "low", f"hard-coded number {v!r} in arithmetic"


TRANSFORMS = {"interp1d", "interp", "interpn", "resample", "savgol_filter",
              "gaussian_filter", "gaussian_filter1d", "convolve", "cumsum",
              "cumulative_trapezoid", "cumtrapz", "trapz", "trapezoid", "gradient",
              "diff", "polyfit", "polyval", "curve_fit", "minimize", "least_squares",
              "fft", "smooth"}
REDUCTIONS = {"sum", "mean", "std", "median", "nanmean", "nansum", "nanstd",
              "average", "var"}
TODO_RE = re.compile(r"\b(TODO|FIXME|XXX|HACK|TBD|WORKAROUND|NOT IMPLEMENTED)\b", re.I)
LOGGER_NAMES = ("logger", "_logger", "logging", "log", "warnings", "print")


@dataclass
class Finding:
    rule: str
    severity: str
    file: str
    line: int
    qualname: str
    message: str
    snippet: str = ""

    @property
    def category(self) -> str:
        return RULES.get(self.rule, "")


class Collector:
    def __init__(self, disabled: set[str]):
        self.items: list[Finding] = []
        self._seen: set[tuple] = set()
        self.disabled = disabled

    def add(self, rule, severity, file, line, qualname, message, snippet=""):
        if rule in self.disabled:
            return
        key = (rule, file, line, message)
        if key in self._seen:
            return
        self._seen.add(key)
        self.items.append(Finding(rule, severity, file, int(line), qualname,
                                  message, snippet.strip()[:160]))


# --------------------------------------------------------- python scanning --
def usrc(node) -> str:
    try:
        return ast.unparse(node)
    except Exception:  # pragma: no cover
        return "<?>"


def is_num(node) -> bool:
    return (isinstance(node, ast.Constant) and isinstance(node.value, (int, float))
            and not isinstance(node.value, bool))


def is_units_key(node) -> bool:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value.endswith("@units") or node.value in ("unit", "units")
    if isinstance(node, ast.JoinedStr) and node.values:
        last = node.values[-1]
        return (isinstance(last, ast.Constant) and isinstance(last.value, str)
                and last.value.endswith("@units"))
    return False


class PyScanner(ast.NodeVisitor):
    def __init__(self, rel: str, lines: list[str], col: Collector, store: dict):
        self.rel, self.lines, self.col, self.store = rel, lines, col, store
        self.prefix = ""
        self.cur = ""
        self.ranges: list[tuple[int, int, str]] = []

    # helpers
    def q(self):
        return self.cur or "<module>"

    def snip(self, node):
        ln = getattr(node, "lineno", 0)
        return self.lines[ln - 1] if 0 < ln <= len(self.lines) else ""

    def add(self, rule, sev, node, msg):
        self.col.add(rule, sev, self.rel, getattr(node, "lineno", 0), self.q(),
                     msg, self.snip(node))

    def unit_value(self, value, node, via):
        for rule, sev, msg in classify_unit(value):
            self.add(rule, sev, node, f"{msg} ({via})")

    def record_key(self, s, node):
        if (isinstance(s, str) and len(s) >= 4 and re.fullmatch(r"[A-Za-z0-9_ \-]+", s)):
            self.store["keys"][re.sub(r"[\s_\-]", "", s).lower()].append(
                (s, self.rel, getattr(node, "lineno", 0)))

    # scopes
    def visit_ClassDef(self, node):
        q = self.prefix + node.name
        self.ranges.append((node.lineno, node.end_lineno, q))
        old = (self.prefix, self.cur)
        self.cur, self.prefix = q, q + "."
        self.generic_visit(node)
        self.prefix, self.cur = old

    def visit_FunctionDef(self, node):
        q = self.prefix + node.name
        self.ranges.append((node.lineno, node.end_lineno, q))
        # numeric defaults
        for d in list(node.args.defaults) + [x for x in node.args.kw_defaults if x]:
            if is_num(d) and d.value not in IGNORE_NUMBERS:
                self.col.add("N03", "info", self.rel, d.lineno, q,
                             f"numeric default parameter {d.value!r} in {node.name}()",
                             self.snip(d))
        body = node.body
        if body and isinstance(body[0], ast.Expr) and isinstance(
                getattr(body[0], "value", None), ast.Constant) and isinstance(
                body[0].value.value, str):
            body = body[1:]
        decorators = " ".join(usrc(d) for d in node.decorator_list)
        if body and all(isinstance(s, ast.Pass) or (
                isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant)
                and s.value.value is Ellipsis) for s in body):
            sev = "info" if "abstract" in decorators else "medium"
            self.col.add("E04", sev, self.rel, node.lineno, q,
                         f"{node.name}() has an empty body (pass / ...)",
                         self.snip(node))
        elif (len(body) == 1 and isinstance(body[0], ast.Raise) and body[0].exc is not None
              and "NotImplementedError" in usrc(body[0].exc)):
            self.col.add("E04", "info", self.rel, node.lineno, q,
                         f"{node.name}() raises NotImplementedError", self.snip(node))
        if len(body) >= 1:  # for K03 / K04
            src = "\n".join(usrc(s) for s in body)
            if src.count("\n") >= 2 and not (node.name.startswith("__") and node.name.endswith("__")):
                self.store["funcs"][node.name].append(
                    (self.rel, q, node.lineno, hashlib.md5(src.encode()).hexdigest(), src))
        old = (self.prefix, self.cur)
        self.cur, self.prefix = q, q + ".<locals>."
        self.generic_visit(node)
        self.prefix, self.cur = old

    visit_AsyncFunctionDef = visit_FunctionDef

    # statements
    def visit_Expr(self, node):
        v = node.value
        if isinstance(v, ast.Constant) and isinstance(v.value, str):  # docstring
            if TODO_RE.search(v.value):
                self.add("E05", "low", node, "TODO/TBD marker in docstring")
            return
        self.generic_visit(node)

    def visit_Assert(self, node):
        self.add("E03", "low", node, "assert is stripped under python -O; "
                                      "do not rely on it for validation")
        self.generic_visit(node)

    def _targets_of(self, node):
        return node.targets if isinstance(node, ast.Assign) else [node.target]

    def visit_Assign(self, node):
        self._assign(node)
        self.generic_visit(node)

    def visit_AnnAssign(self, node):
        if node.value is not None:
            self._assign(node)
        self.generic_visit(node)

    def _assign(self, node):
        value = node.value
        for t in self._targets_of(node):
            if isinstance(t, ast.Name) and re.search(r"unit", t.id, re.I) and isinstance(value, ast.Dict):
                for v in value.values:
                    if isinstance(v, ast.Constant) and isinstance(v.value, str):
                        self.unit_value(v.value, v, f"value in {t.id}")
            if isinstance(t, ast.Subscript) and is_units_key(t.slice):
                if isinstance(value, ast.Constant) and isinstance(value.value, str):
                    self.add("U06", "medium", node,
                             f"unit label {value.value!r} is hard-coded, not read from the file")
                    self.unit_value(value.value, node, "assigned to a @units key")

    def visit_Try(self, node):
        for h in node.handlers:
            tname = usrc(h.type) if h.type is not None else None
            body = h.body
            swallowed = all(
                isinstance(s, (ast.Pass, ast.Continue, ast.Break))
                or (isinstance(s, ast.Return) and (s.value is None or isinstance(s.value, ast.Constant)))
                or (isinstance(s, ast.Expr) and isinstance(s.value, ast.Call)
                    and usrc(s.value.func).split(".")[0] in LOGGER_NAMES)
                for s in body)
            if h.type is None:
                self.add("E01", "medium", h, "bare 'except:' hides every error type")
            elif swallowed:
                broad = tname in ("Exception", "BaseException")
                self.add("E02", "medium" if broad else "low", h,
                         f"except {tname}: error is swallowed and execution continues")
            if h.type is None and swallowed:
                self.add("E02", "medium", h, "exception swallowed in bare except")
        self.generic_visit(node)

    visit_TryStar = visit_Try

    # expressions
    def visit_BinOp(self, node):
        if isinstance(node.op, (ast.Mult, ast.Div, ast.FloorDiv, ast.Pow)):
            for side in (node.left, node.right):
                if is_num(side) and side.value not in IGNORE_NUMBERS:
                    rule, sev, msg = number_hint(side.value)
                    self.add(rule, sev, node, msg)
        self.generic_visit(node)

    def visit_BoolOp(self, node):
        if isinstance(node.op, ast.Or) and node.values:
            last = node.values[-1]
            if isinstance(last, ast.Constant) and last.value not in (None, "", 0, False) \
                    and not isinstance(last.value, bool):
                self.add("D03", "low", node, f"'... or {last.value!r}' silently replaces "
                                              "a falsy value (incl. 0 / empty)")
        self.generic_visit(node)

    def visit_Compare(self, node):
        for op, comp in zip(node.ops, node.comparators):
            consts = []
            if isinstance(op, (ast.Eq, ast.NotEq)) and isinstance(comp, ast.Constant):
                consts = [comp.value]
            elif isinstance(op, (ast.In, ast.NotIn)) and isinstance(comp, (ast.List, ast.Tuple, ast.Set)):
                consts = [e.value for e in comp.elts if isinstance(e, ast.Constant)]
            for c in consts:
                if isinstance(c, str) and len(c) >= 3 and re.search(r"[A-Za-z]", c) \
                        and not isinstance(node.left, ast.Constant) \
                        and "__main__" not in c:
                    self.add("D05", "low", node,
                             f"exact, case-sensitive match against {c!r}; other spellings fall through")
        self.generic_visit(node)

    def visit_Attribute(self, node):
        if usrc(node).endswith("timezone.utc"):
            self.add("T01", "medium", node, "UTC timezone assumed/attached here")
        self.generic_visit(node)

    def visit_Constant(self, node):
        v = node.value
        if isinstance(v, str) and ("ENTRY[" in v or v.startswith("/ENTRY")):
            self.add("C01", "info", node, f"NeXus path literal {v[:60]!r} in code "
                                          "(mapping outside the config files)")

    def visit_Dict(self, node):
        seen = {}
        for k, v in zip(node.keys, node.values):
            if isinstance(k, ast.Constant):
                key = (type(k.value), k.value)
                if key in seen:
                    self.add("K01", "high", k, f"duplicate dict key {k.value!r} "
                             f"(first at line {seen[key]}); the last value silently wins")
                else:
                    seen[key] = k.lineno
                if isinstance(k.value, str):
                    self.record_key(k.value, k)
            if k is not None and is_units_key(k) and isinstance(v, ast.Constant) \
                    and isinstance(v.value, str):
                self.unit_value(v.value, v, "dict entry for a unit key")
        self.generic_visit(node)

    def visit_Subscript(self, node):
        sl = node.slice
        if isinstance(sl, ast.Slice) and isinstance(sl.step, ast.UnaryOp) \
                and isinstance(sl.step.op, ast.USub) and is_num(sl.step.operand) \
                and sl.step.operand.value == 1:
            self.add("A01", "medium", node, "reversed with [::-1]: check the paired axis/signal")
        if isinstance(sl, ast.Constant) and isinstance(sl.value, str):
            self.record_key(sl.value, node)
        if (isinstance(sl, ast.Constant) and sl.value == 0 and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Attribute) and node.value.func.attr == "split"
                and node.value.args and isinstance(node.value.args[0], ast.Constant)
                and node.value.args[0].value == "."):
            self.add("N05", "low", node, "split('.')[0] drops the fractional part")
        self.generic_visit(node)

    def visit_Call(self, node):
        name = usrc(node.func)
        last = name.split(".")[-1]

        # rounding
        if last in ("round", "around"):
            nd = node.args[1] if len(node.args) > 1 else next(
                (k.value for k in node.keywords if k.arg in ("ndigits", "decimals")), None)
            if isinstance(nd, ast.Constant) and isinstance(nd.value, int):
                sev = "medium" if nd.value <= 3 else "low"
                self.add("N04", sev, node, f"rounded to {nd.value} decimals: "
                                           "loses precision for fine steps / small values")
        # flips
        if last in ("flip", "flipud", "fliplr") or (last == "reverse" and not node.args):
            self.add("A01", "medium", node, f"{name}(): axis reversed; "
                                           "verify the paired signal is reordered too")
        # defaults
        if isinstance(node.func, ast.Attribute) and node.func.attr == "get" and node.args:
            a0 = node.args[0]
            if isinstance(a0, ast.Constant) and isinstance(a0.value, str):
                self.record_key(a0.value, node)
            if len(node.args) == 2 and isinstance(node.args[1], ast.Constant) \
                    and node.args[1].value is not None:
                dflt = node.args[1].value
                target = usrc(node.func.value)
                if re.search(r"kwargs|options|params|settings", target):
                    self.add("D02", "medium", node, f"option default {dflt!r} for "
                             f"{usrc(a0)} read from {target}: check docs match this default")
                else:
                    sev = "info" if dflt in ("", 0, False) else "low"
                    self.add("D01", sev, node, f".get({usrc(a0)}, {dflt!r}): silent fallback value")
        if name == "getattr" and len(node.args) == 3 and isinstance(node.args[2], ast.Constant) \
                and node.args[2].value is not None:
            self.add("D01", "low", node, f"getattr default {node.args[2].value!r}")
        # fabricated arrays
        if last in ("ones", "ones_like", "full", "full_like"):
            self.add("D04", "medium", node, f"{name}(): constant array; if it stands in for "
                                           "missing data it fabricates values")
        elif last in ("zeros", "zeros_like"):
            self.add("D04", "low", node, f"{name}(): zero-filled array")
        # numerical processing
        if last in TRANSFORMS:
            self.add("A02", "low", node, f"{name}(): data-changing operation; "
                                         "is the processing recorded in the output?")
        elif last in REDUCTIONS:
            self.add("A03", "info", node, f"{name}(): reduction; check mean vs sum "
                                          "and the naming in the NeXus output")
        # casts / truncations
        if name == "int" and node.args and isinstance(node.args[0], ast.Call) \
                and usrc(node.args[0].func) == "float":
            self.add("N05", "low", node, "int(float(x)) truncates toward zero")
        if name == "min" and len(node.args) == 2 and all(
                isinstance(a, ast.Call) and usrc(a.func) == "len" for a in node.args):
            self.add("N06", "low", node, "min(len(a), len(b)): arrays of different "
                                         "length are silently truncated")
        # time handling
        if last == "datetime" and len(node.args) >= 3 and all(
                is_num(a) for a in node.args[:3]) and [a.value for a in node.args[:3]] == [1, 1, 1]:
            self.add("T02", "medium", node, "sentinel date 0001-01-01 used as fallback")
        if last in ("timezone", "tzoffset", "tzutc") or any(k.arg == "tzinfo" for k in node.keywords):
            self.add("T01", "medium", node, f"timezone constructed/attached via {name}()")
        if last in ("utcnow", "utcfromtimestamp", "fromtimestamp", "strptime", "now"):
            self.add("T03", "info", node, f"{name}(): naive or format-dependent time handling")
        # units passed as keywords
        for k in node.keywords:
            if k.arg in ("unit", "units") and isinstance(k.value, ast.Constant) \
                    and isinstance(k.value.value, str):
                self.unit_value(k.value.value, node, f"keyword {k.arg}=")
        if name == "print":
            self.add("M01", "info", node, "print() in library code")
        self.generic_visit(node)


def scan_comments(rel, text, ranges, col):
    try:
        toks = tokenize.generate_tokens(io.StringIO(text).readline)
        for tok in toks:
            if tok.type == tokenize.COMMENT and TODO_RE.search(tok.string):
                ln = tok.start[0]
                inner = [r for r in ranges if r[0] <= ln <= r[1]]
                qn = min(inner, key=lambda r: r[1] - r[0])[2] if inner else "<module>"
                col.add("E05", "low", rel, ln, qn, tok.string.strip()[:100], tok.line)
    except (tokenize.TokenError, IndentationError):
        pass


# ----------------------------------------------------------- config files ---
JSON_LINE = re.compile(r'^\s*"(?P<key>(?:[^"\\]|\\.)*)"\s*:\s*(?P<val>.*?)\s*,?\s*$')
YAML_LINE = re.compile(r'^\s*-?\s*(?P<key>[A-Za-z0-9_@/\[\]\. \-]+?)\s*:\s+(?P<val>[^#]+?)\s*(#.*)?$')


def scan_json(rel, text, col):
    dups: list[str] = []

    def hook(pairs):
        c = Counter(k for k, _ in pairs)
        dups.extend(k for k, n in c.items() if n > 1)
        return dict(pairs)

    try:
        json.loads(text, object_pairs_hook=hook)
    except json.JSONDecodeError as e:
        col.add("C05", "medium", rel, e.lineno, "<config>", f"JSON parse error: {e.msg}")
        return
    lines = text.splitlines()
    for key in set(dups):
        hits = [i for i, l in enumerate(lines, 1) if f'"{key}"' in l and ":" in l]
        for i in hits[1:]:
            col.add("C04", "high", rel, i, "<config>",
                    f"duplicate JSON key {key[:80]!r}: the last occurrence wins", lines[i - 1])
    for i, line in enumerate(lines, 1):
        m = JSON_LINE.match(line)
        if not m:
            continue
        key, raw = m["key"], m["val"]
        if raw.startswith(("{", "[")):
            continue
        try:
            v = json.loads(raw)
        except ValueError:
            continue
        unit_key = key.endswith("@units") or key in ("unit", "units")
        if isinstance(v, str):
            if v.startswith(("@", "!@")):
                continue  # lookup into parsed data / ELN: mapping, not hard-coded
            if unit_key:
                col.add("U07", "low", rel, i, "<config>",
                        f"unit {v!r} is fixed for every file (key {key[-50:]!r})", line)
                for rule, sev, msg in classify_unit(v):
                    col.add(rule, sev, rel, i, "<config>", msg + " (config)", line)
            else:
                col.add("C02", "info", rel, i, "<config>",
                        f"hard-coded value {v[:50]!r} for {key[-60:]!r}", line)
        elif isinstance(v, (int, float)) and not isinstance(v, bool):
            col.add("C03", "low", rel, i, "<config>",
                    f"hard-coded number {v!r} for {key[-60:]!r}", line)


def scan_yaml(rel, text, col):
    for i, line in enumerate(text.splitlines(), 1):
        m = YAML_LINE.match(line)
        if m and m["key"].strip() in ("unit", "units"):
            val = m["val"].strip().strip("'\"")
            if val.startswith(("@", "!@", "{", "[")):
                continue
            for rule, sev, msg in classify_unit(val):
                col.add(rule, sev, rel, i, "<config>", msg + " (yaml)", line)


# ------------------------------------------------------------ cross-file ----
def cross_file(store, col):
    for norm, uses in store["keys"].items():
        variants = defaultdict(list)
        for s, f, ln in uses:
            variants[s].append((f, ln))
        if len(variants) < 2:
            continue
        ordered = sorted(variants.items(), key=lambda kv: -len(kv[1]))
        desc = "; ".join(f"{s!r} x{len(l)} (e.g. {l[0][0]}:{l[0][1]})" for s, l in ordered[:5])
        minor = ordered[-1][1][0]
        col.add("K02", "low", minor[0], minor[1], "<keys>",
                f"same key spelled differently: {desc}")
    ignore = {"_parse", "parse", "matches_file", "detect_version", "_build_parsed_spectra"}
    for fname, defs in store["funcs"].items():
        if len(defs) < 2 or fname in ignore:
            continue
        hashes = {d[3] for d in defs}
        first = defs[0]
        locs = ", ".join(f"{d[0]}:{d[2]}" for d in defs[:6])
        if len(hashes) == 1:
            col.add("K04", "info", first[0], first[2], first[1],
                    f"{fname}() duplicated verbatim in: {locs}")
        else:
            best = max(difflib.SequenceMatcher(None, first[4], d[4]).ratio()
                       for d in defs[1:])
            if best >= 0.6:
                col.add("K03", "low", first[0], first[2], first[1],
                        f"{fname}() exists in similar but different copies "
                        f"(similarity up to {best:.2f}): {locs}. Compare behaviour")


# -------------------------------------------------------------- driver ------
def iter_files(root: Path, include_tests: bool, include_nomad: bool):
    base = root / "src" if (root / "src").is_dir() else root
    roots = [base] + ([root / "tests"] if include_tests and (root / "tests").is_dir() else [])
    for b in roots:
        for f in sorted(b.rglob("*")):
            if not f.is_file():
                continue
            rel = f.relative_to(root)
            parts = rel.parts
            if any(p in SKIP_DIRS or p.endswith(".egg-info") for p in parts):
                continue
            if not include_nomad and "nomad" in parts:
                continue
            if f.suffix in (".py", ".json", ".yaml", ".yml"):
                yield f, str(rel).replace("\\", "/")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("root", nargs="?", default=".")
    ap.add_argument("--out", default="unorthodox_findings.csv")
    ap.add_argument("--inventory", help="path to inventory.csv to annotate")
    ap.add_argument("--min-severity", choices=list(SEV), default="info")
    ap.add_argument("--disable", default="", help="comma-separated rule ids to skip")
    ap.add_argument("--include-tests", action="store_true")
    ap.add_argument("--include-nomad", action="store_true",
                    help="also scan nomad/ example uploads (skipped by default)")
    ap.add_argument("--list-rules", action="store_true")
    a = ap.parse_args()

    if a.list_rules:
        for k, v in RULES.items():
            print(f"{k}  {v}")
        return
    root = Path(a.root).resolve()
    col = Collector({r.strip() for r in a.disable.split(",") if r.strip()})
    store = {"keys": defaultdict(list), "funcs": defaultdict(list)}
    n_files = Counter()

    for f, rel in iter_files(root, a.include_tests, a.include_nomad):
        text = f.read_text(errors="ignore")
        n_files[f.suffix] += 1
        if f.suffix == ".py":
            try:
                tree = ast.parse(text)
            except SyntaxError as e:
                col.add("C05", "medium", rel, e.lineno or 1, "<module>", f"syntax error: {e.msg}")
                continue
            sc = PyScanner(rel, text.splitlines(), col, store)
            sc.visit(tree)
            scan_comments(rel, text, sc.ranges, col)
        elif f.suffix == ".json":
            scan_json(rel, text, col)
        else:
            scan_yaml(rel, text, col)
    cross_file(store, col)

    items = [x for x in col.items if SEV[x.severity] >= SEV[a.min_severity]]
    items.sort(key=lambda x: (-SEV[x.severity], x.file, x.line))
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rule", "severity", "category", "file", "line", "qualname",
                    "message", "snippet", "reviewed", "notes"])
        for x in items:
            w.writerow([x.rule, x.severity, x.category, x.file, x.line, x.qualname,
                        x.message, x.snippet, "", ""])

    print(f"scanned: {dict(n_files)}   findings: {len(items)}  -> {a.out}")
    by_sev = Counter(x.severity for x in items)
    print("by severity:", {s: by_sev.get(s, 0) for s in ("high", "medium", "low", "info")})
    print("by rule:")
    for r, n in Counter(x.rule for x in items).most_common():
        print(f"  {r}  {n:4d}  {RULES.get(r, '')}")
    print("files with most medium/high findings:")
    hot = Counter(x.file for x in items if SEV[x.severity] >= 2)
    for f, n in hot.most_common(15):
        print(f"  {n:4d}  {f}")

    if a.inventory:
        with open(a.inventory, newline="") as fh:
            inv = list(csv.DictReader(fh))
        by_file = defaultdict(list)
        for i, row in enumerate(inv):
            by_file[row["file"]].append(i)
        stats = defaultdict(lambda: {"n": 0, "max": -1, "rules": Counter()})
        for x in items:
            best, span = None, 10 ** 9
            for i in by_file.get(x.file, []):
                row = inv[i]
                try:
                    s, e = int(float(row["line_start"])), int(float(row["line_end"]))
                except (ValueError, TypeError):
                    continue
                if s <= x.line <= e and (e - s) < span and row["kind"] != "module":
                    best, span = i, e - s
            if best is None:  # fall back to the module row
                for i in by_file.get(x.file, []):
                    if inv[i]["kind"] in ("module", "config"):
                        best = i
            if best is not None:
                st = stats[best]
                st["n"] += 1
                st["max"] = max(st["max"], SEV[x.severity])
                st["rules"][x.rule] += 1
        names = {v: k for k, v in SEV.items()}
        out = a.inventory.replace(".csv", "") + "_with_findings.csv"
        cols = list(inv[0].keys()) + ["scan_findings", "scan_max_severity", "scan_rules"]
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for i, row in enumerate(inv):
                st = stats.get(i)
                row["scan_findings"] = st["n"] if st else 0
                row["scan_max_severity"] = names[st["max"]] if st else ""
                row["scan_rules"] = ";".join(f"{r}x{n}" for r, n in sorted(st["rules"].items())) if st else ""
                w.writerow(row)
        print(f"annotated ledger -> {out}")


if __name__ == "__main__":
    sys.exit(main())
