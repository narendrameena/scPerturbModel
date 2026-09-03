#!/usr/bin/env python3
"""Verify every headline number in RESULTS.md and MANUSCRIPT.md against source.

This project has been corrected repeatedly -- the interaction estimator, the
matched null, the leave-one-context-out prior, the plate-14 exclusion, the
withdrawn mechanism claim, the reservoir-sampled bootstrap. Each correction
re-ran some analyses and not others, so the live risk is a number in the prose
that was computed under a superseded code path and is now simply wrong. Prose
does not recompute itself.

This script recomputes each headline quantity from the table the analysis wrote,
and checks the value actually appears in the documents. It is deliberately dumb:
it does not re-derive anything, it only asks whether what we claim matches what
the current pipeline produced. A FAIL means either the prose is stale or the
table is, and both need looking at.

Run it before any submission, and after any change to an estimator.

Outputs: results/tables/number_audit.csv (exit status is non-zero on failure)
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
LOGS = ROOT / "logs"
DOCS = [ROOT / "RESULTS.md", ROOT / "MANUSCRIPT.md"]


def newest_log(pattern, regex, group=1):
    """Pull a value out of the most recent matching SLURM log."""
    for f in sorted(LOGS.glob(pattern), key=lambda p: p.stat().st_mtime,
                    reverse=True):
        m = re.search(regex, f.read_text(errors="ignore"))
        if m:
            return float(m.group(group)), f.name
    return None, None


def table(name):
    f = TAB / name
    return pd.read_csv(f) if f.exists() else None


def check(rows, claim, computed, quoted_as, docs_text, tol=0.02, unit=""):
    """Record whether `quoted_as` appears in the prose and matches `computed`."""
    present = quoted_as in docs_text
    ok = present
    note = "" if present else f"'{quoted_as}' not found in prose"
    rows.append({"claim": claim, "computed": computed, "quoted": quoted_as,
                 "in_prose": present, "status": "OK" if ok else "CHECK",
                 "note": note})


def main():
    docs = "\n".join(p.read_text() for p in DOCS if p.exists())
    rows = []

    # --- Tahoe replicate structure and interaction share -------------------
    t = table("tahoe_replicate_structure.csv")
    if t is not None:
        v = float(t[t["filter"] == "pass_filter == full"].frac_replicated.iloc[0])
        # two decimals, so the check cannot be satisfied by a rounded form of
        # the same number appearing elsewhere in the prose -- which is how a
        # blanket string replacement once corrupted 13.52% to 13.46% while this
        # audit still passed
        check(rows, "Tahoe replication (atlas QC filter)", round(v, 4),
              f"{v*100:.2f}%", docs)
    t = table("tahoe_true_replicates.csv")
    if t is not None:
        tr = t[t.pairing.str.startswith("true replicate")]
        if len(tr):
            s = float(tr.share.iloc[0])
            check(rows, "Tahoe true-replicate interaction share", round(s, 4),
                  f"{s*100:.1f}%", docs)
            check(rows, "Tahoe true-replicate CI low", round(float(tr.share_lo.iloc[0]), 4),
                  f"{tr.share_lo.iloc[0]*100:.1f}", docs)
        xd = t[t.pairing.str.startswith("cross-dose")]
        if len(xd):
            s = float(xd.share.iloc[0])
            check(rows, "Tahoe cross-dose share", round(s, 4),
                  f"{s*100:.1f}%", docs)

    # --- decomposition shares from the recorded logs ----------------------
    for lab, pat, claim in (
            ("LINCS p1", "lincs_p1_*.out", "LINCS phase 1 interaction share"),
            ("LINCS p2", "lincs_mech_*.out", "LINCS phase 2 interaction share")):
        v, src = newest_log(pat, r"context x perturbation\s+[\d.]+\s+\((\d+)%")
        if v is not None:
            check(rows, claim, v, f"{int(v)}%", docs)
    v, src = newest_log("prism_*.out",
                        r"interaction ([\d.]+) \((\d+)%\)", group=2)
    if v is not None:
        check(rows, "PRISM interaction share", v, f"{int(v)}", docs)

    # --- cross-laboratory quantities --------------------------------------
    s = table("cross_lab_summary.csv")
    if s is not None:
        d = dict(zip(s.quantity, s.value))
        for key, claim, fmt in (
                ("PRISM replicate ceiling", "PRISM within-lab ceiling", "{:.3f}"),
                ("GDSC1 vs GDSC2 ceiling", "GDSC within-lab ceiling", "{:.3f}"),
                ("cross-lab PRISM vs GDSC (all ID-matched)",
                 "cross-lab median r", "{:.3f}"),
                ("reproducible fraction (all)", "reproducible fraction", "{:.1%}"),
                ("share of loss from changing laboratory",
                 "laboratory share of loss", "{:.0%}")):
            if key in d and np.isfinite(d[key]):
                check(rows, claim, round(float(d[key]), 4),
                      fmt.format(d[key]), docs)

    # --- transcription cross-lab arm (added after a stale value slipped in
    # when the sciPlex3 rows were pooled into the headline by accident) ------
    tx = table("cross_lab_transcription.csv")
    if tx is not None and len(tx):
        wl = tx[tx.comparison.str.contains("within-lab")
                & ~tx.comparison.str.contains("validated")].r.median()
        xl = tx[tx.comparison.str.contains("Tahoe vs LINCS")
                & ~tx.comparison.str.contains("validated")].r.median()
        # The transcriptional arm has no cross-laboratory comparison at all:
        # after identity validation no (line, compound) pair is shared between
        # Tahoe and LINCS. Checking prose for a value that does not exist would
        # report a permanent failure and train the reader to ignore this audit.
        if np.isfinite(wl) and wl > 0 and np.isfinite(xl):
            check(rows, "transcription reproducible fraction",
                  round(float(xl / wl), 3), f"{xl/wl*100:.0f}%", docs)

    # --- predictor blocks --------------------------------------------------
    e = table("expression_architecture.csv")
    if e is not None:
        for col, claim in (("r2_expression", "expression CV R2"),
                           ("r2_protein", "protein CV R2"),
                           ("r2_lineage", "lineage CV R2"),
                           ("r2_cnv", "copy-number CV R2"),
                           ("r2_nonsyn", "nonsynonymous CV R2")):
            if col in e.columns and e[col].notna().any():
                v = float(e[col].median())
                check(rows, claim, round(v, 4), f"{v:+.4f}".rstrip("0"), docs)

    # --- power curve -------------------------------------------------------
    p = table("prism_genotype_power.csv")
    if p is not None:
        for n in (47, 400, 600):
            r = p[p.n_contexts == n]
            if len(r):
                v = float(r.detection_rate.iloc[0])
                check(rows, f"genotype power at {n} contexts", round(v, 3),
                      f"{v*100:.0f}%", docs)

    # --- baseline scaling / benchmark handicap -----------------------------
    a = table("additive_vs_n_lines.csv")
    if a is not None:
        g = a.groupby("k").r_de100.mean()
        for n in (2, 18):
            if n in g.index:
                check(rows, f"baseline r at {n} contexts", round(float(g[n]), 4),
                      f"{g[n]:.3f}", docs)
        if 2 in g.index and 45 in g.index:
            h = (g[45] - g[2]) / g[2] * 100
            check(rows, "handicap, 2 vs 45 contexts", round(float(h), 1),
                  f"{h:.0f}%", docs)

    # --- identity stringency ------------------------------------------------
    st = table("identity_stringency.csv")
    if st is not None and len(st):
        top5 = st[st.criterion.str.contains("top 5", na=False)]
        if len(top5):
            check(rows, "identity 'top 5%' lines retained",
                  int(top5.n_lines.iloc[0]), str(int(top5.n_lines.iloc[0])),
                  docs)

    # --- matching ladder ----------------------------------------------------
    m = table("matching_strategy.csv")
    if m is not None and len(m):
        g = m.groupby("strategy").r.median()
        for key, claim in ((("identifier"), "identifier matching r"),
                           (("best hit"), "best-hit r")):
            hit = [k for k in g.index if key in k]
            if hit:
                v = float(g[hit[0]])
                check(rows, claim, round(v, 3), f"{v:.3f}", docs)

    # --- simulation ---------------------------------------------------------
    sim = table("estimator_simulation.csv")
    if sim is not None and len(sim):
        z = sim[(sim.sweep == "share") & (sim.true_share == 0.0)]
        for est, claim in (("residual variance", "residual variance at true 0"),
                           ("pooled batch", "pooled batch at true 0")):
            s_ = z[z.estimator == est].estimate
            if len(s_):
                v = float(s_.mean())
                check(rows, claim, round(v, 3), f"{v*100:.0f}%", docs)

    A = pd.DataFrame(rows)
    A.to_csv(TAB / "number_audit.csv", index=False)
    n_ok = int((A.status == "OK").sum())
    print(f"{n_ok}/{len(A)} headline numbers verified present and current\n")
    print(A[["claim", "computed", "quoted", "status"]].to_string(index=False))
    bad = A[A.status != "OK"]
    if len(bad):
        print(f"\n{len(bad)} need attention — either the prose is stale or the "
              f"table is:")
        for r in bad.itertuples():
            print(f"  {r.claim}: computed {r.computed}, {r.note}")
        return 1
    print("\nAll checked numbers match the current tables.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
