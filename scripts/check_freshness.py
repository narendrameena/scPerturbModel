#!/usr/bin/env python3
"""Which results are older than the code or data that produced them.

Every wrong number found in the 2026-09-07/08 audits was a *stale* number: a
table generated before an input was regenerated, or before the script that writes
it was corrected. Two cases cost real work.

  `three_platform_mechanism_cdi.csv` was written 2026-08-31; three of its four
  inputs were regenerated on 09-01 and 09-03 when the corrected estimator landed.
  Nothing re-ran it, so §16 and §31 quoted values no saved artefact reproduced.

  `cross_lab_transcription.py` was corrected on 09-03 in a way that silently
  broke its keying. Its outputs were left at 09-01, so for five days the repo
  held a table describing neither the old code nor the new.

Neither is detectable by reading the output — both were plausible. Both are
trivially detectable by comparing timestamps, which is what this does:

    output is STALE if it is older than
      * the script that writes it,                    <- catches the second case
      * the specific library FUNCTIONS it imports,    <- catches estimator changes
      * or any table it reads.                        <- catches the first case

Library dependencies resolve to functions, not files. `delta_eval.py` gained a
new function on 09-03 without touching the two that most scripts call; at file
granularity that marked nine outputs stale which were not. `git log -L :fn:file`
gives the last change to each definition, so only callers of a changed function
are flagged.

Dependencies are read from the source, not declared by hand, so a new script is
covered the moment it is written.

    python scripts/check_freshness.py            # report
    python scripts/check_freshness.py --rerun    # print the commands to fix it
    python scripts/check_freshness.py --quiet    # exit 1 if anything is stale

`tests/test_celldrug.py` runs this over the tables backing manuscript numbers, so
a stale headline fails the suite rather than reaching a reviewer.

KNOWN FALSE POSITIVE. A code file's time is the later of its filesystem mtime and
its last commit, because a checkout resets mtimes and would otherwise hide a real
change (`cross_lab_reproducibility.py` carries an mtime an hour before the commit
that wrote it). The cost is that the ordinary workflow -- run the script, then
commit it -- leaves the commit a minute or two after the output, and that reads as
stale. Gaps of minutes are usually this; gaps of days are usually real. Results
are sorted by gap size for that reason.

The durable fix is to record the hash of the code that produced each output
rather than infer it from timestamps. That is not implemented.
"""
from __future__ import annotations

import argparse
import ast
import re
import sys
import subprocess
from datetime import datetime
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
FIGROOT = ROOT / "results" / "figures"
SRC = ROOT / "src" / "perturbmodel"

# Artefacts written in one call land microseconds apart; save_figure writes the
# PNG then its source CSV. 60 s is far below any real staleness (the two cases
# above were days) and far above same-run jitter (max observed 37 s).
TOLERANCE_S = 60


def mtime(p: Path) -> float:
    return p.stat().st_mtime


def _strings(node) -> list[str]:
    """String literals in an expression, including f-string fixed parts."""
    out = []
    for n in ast.walk(node):
        if isinstance(n, ast.Constant) and isinstance(n.value, str):
            out.append(n.value)
    return out


@lru_cache(maxsize=None)
def commit_mtime(rel: str) -> float | None:
    """When git last committed this file."""
    try:
        r = subprocess.run(["git", "log", "-1", "--format=%ct", "--", rel],
                           cwd=ROOT, capture_output=True, text=True, timeout=30)
        t = r.stdout.strip()
        return float(t) if t.isdigit() else None
    except Exception:
        return None


def code_mtime(p: Path) -> float:
    """A code file's true last-change time.

    The later of its filesystem mtime and its last commit. Neither alone is
    right: a checkout resets mtimes to the checkout time, so a file committed at
    12:30 can carry an mtime of 11:28 (`cross_lab_reproducibility.py` does), and
    conversely an uncommitted edit has no commit time at all.
    """
    t = mtime(p)
    c = commit_mtime(str(p.relative_to(ROOT)))
    return max(t, c) if c else t


@lru_cache(maxsize=None)
def symbol_mtime(rel: str, name: str) -> float | None:
    """When this definition last changed, from git. None if not determinable."""
    try:
        r = subprocess.run(
            ["git", "log", "-1", "-L", f":{name}:{rel}", "--format=%ct"],
            cwd=ROOT, capture_output=True, text=True, timeout=30)
        if r.returncode != 0:
            return None
        for line in r.stdout.splitlines():
            if line.strip().isdigit():
                return float(line.strip())
    except Exception:
        return None
    return None


def parse_script(path: Path) -> dict:
    """Outputs, table inputs and imported project symbols, read from the source."""
    try:
        tree = ast.parse(path.read_text())
    except SyntaxError:
        return {"writes": set(), "reads": set(), "modules": set(),
                "figures": set()}
    writes, reads, modules, figures = set(), set(), set(), set()
    for n in ast.walk(tree):
        if isinstance(n, (ast.Import, ast.ImportFrom)):
            mod = getattr(n, "module", None) or ""
            names = [a.name for a in n.names]
            if mod.startswith("perturbmodel"):
                # (module, symbol) so staleness resolves per definition;
                # symbol None means the whole module is a dependency
                for a in names:
                    modules.add((mod, a))
            for a in names:
                if a.startswith("perturbmodel"):
                    modules.add((a, None))
        if not isinstance(n, ast.Call):
            continue
        fn = n.func
        name = fn.attr if isinstance(fn, ast.Attribute) else \
            (fn.id if isinstance(fn, ast.Name) else "")
        if name == "to_csv":
            writes |= {s for s in _strings(n) if s.endswith((".csv", ".csv.gz"))}
        elif name == "read_csv":
            reads |= {s for s in _strings(n) if s.endswith((".csv", ".csv.gz"))}
        elif name == "save_figure" and len(n.args) >= 2:
            # save_figure(fig, name, out_root) -> out_root/<name>/<name>.png
            for lit in _strings(n.args[1]):
                if lit and "/" not in lit:
                    figures.add(lit)
    return {"writes": writes, "reads": reads, "modules": modules,
            "figures": figures}


def module_deps(modules: set) -> list[tuple[Path, str | None, float]]:
    """(file, symbol, effective mtime) for each imported project symbol.

    A symbol resolves to when its definition last changed; a bare module import,
    or a symbol git cannot locate (a re-export, a class), falls back to the
    file's mtime, which is conservative.
    """
    out = []
    for mod, sym in modules:
        rel = Path(*mod.split(".")[1:])
        f = None
        for cand in (SRC / rel.with_suffix(".py"), SRC / rel / "__init__.py"):
            if cand.exists():
                f = cand
                break
        if f is None:
            continue
        t = None
        if sym:
            t = symbol_mtime(str(f.relative_to(ROOT)), sym)
        out.append((f, sym, t if t is not None else code_mtime(f)))
    return out


def resolve_figures(name: str) -> list[Path]:
    """A save_figure name to the PNG it writes, wherever the bundle landed.

    `manuscript_figures.py` calls `save_figure(fig, f"fig{i}", ...)`, whose only
    string literal is "fig". Matching that exactly resolved to a directory that
    does not exist, so fig1..fig10 -- every figure the design paper cites -- were
    invisible to this checker while it reported the repository clean. Treat a
    name that matches no bundle exactly as a prefix instead.
    """
    exact = sorted(FIGROOT.rglob(f"{name}/{name}.png"))
    if exact:
        return exact
    return sorted(p for p in FIGROOT.rglob("*/*.png")
                  if p.parent.name.startswith(name) and p.stem == p.parent.name)


def resolve(pattern: str) -> list[Path]:
    """A literal table name, or the fixed part of an f-string name, to real files."""
    base = Path(pattern).name
    if "{" not in base:
        p = TAB / base
        return [p] if p.exists() else []
    stem = base.split("{")[0]
    return sorted(TAB.glob(f"{stem}*")) if stem else []


def audit() -> list[dict]:
    rows = []
    for script in sorted((ROOT / "scripts").glob("*.py")):
        if script.name == Path(__file__).name:
            continue
        info = parse_script(script)
        if not (info["writes"] or info["figures"]):
            continue
        deps = [(script, None, code_mtime(script))]
        deps += module_deps(info["modules"])
        read_paths = set()
        for r in info["reads"]:
            for p in resolve(r):
                read_paths.add(p)
                deps.append((p, None, mtime(p)))
        outputs = [o for w in info["writes"] for o in resolve(w)]
        outputs += [o for f_ in info["figures"] for o in resolve_figures(f_)]
        for out in outputs:
            if True:
                if out in read_paths:    # --replot: reads what it writes
                    continue
                t = mtime(out)
                newer = [(d, s_, dt) for d, s_, dt in deps
                         if dt > t + TOLERANCE_S and d != out]
                if newer:
                    why = max(newer, key=lambda x: x[2])
                    rows.append({"output": out, "script": script, "out_t": t,
                                 "cause": why[0], "symbol": why[1],
                                 "cause_t": why[2], "n_causes": len(newer)})
    return sorted(rows, key=lambda r: r["cause_t"] - r["out_t"], reverse=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rerun", action="store_true",
                    help="print the commands that would refresh what is stale")
    ap.add_argument("--quiet", action="store_true",
                    help="print nothing on success; exit 1 if anything is stale")
    ap.add_argument("--critical", action="store_true",
                    help="only tables mirrored in docs/source_data/, i.e. those "
                         "backing a number quoted in a manuscript")
    a = ap.parse_args()
    rows = audit()
    if a.critical:
        crit = {p.name for p in (ROOT / "docs" / "source_data").glob("*.csv")}
        rows = [r for r in rows if r["output"].name in crit]
    if not rows:
        if not a.quiet:
            scope = ("Every table backing a manuscript number is"
                     if a.critical else "All results are")
            print(f"{scope} newer than the code and data that produced it.")
        return 0
    if a.rerun:
        for s in sorted({r["script"].name for r in rows}):
            print(f"python scripts/{s}")
        return 1
    if a.quiet:
        return 1
    fmt = "%m-%d %H:%M"
    print(f"{len(rows)} STALE RESULT(S) — older than what produced them:\n")
    for r in rows:
        age = (r["cause_t"] - r["out_t"]) / 86400
        print(f"  {r['output'].name}")
        what = f"{r['cause'].relative_to(ROOT)}" + (
            f"::{r['symbol']}" if r.get("symbol") else "")
        print(f"    written {datetime.fromtimestamp(r['out_t']):{fmt}}  "
              f"but {what} changed "
              f"{datetime.fromtimestamp(r['cause_t']):{fmt}}  "
              f"({age:.1f} d newer"
              + (f", {r['n_causes']} deps newer)" if r["n_causes"] > 1 else ")"))
        print(f"    fix: python scripts/{r['script'].name}")
    print("\n  --rerun prints just the commands.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
