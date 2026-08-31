#!/usr/bin/env python3
"""Run each rejected methodological alternative against the chosen one.

A rationale that argues from principle is an opinion; a rationale that shows the
rejected alternative producing a different number on the same data is evidence.
This script does the latter for every decision recorded in
docs/methodology_rationale.md that is computable, so each entry cites a measured
quantity with a script behind it rather than a claim.

Each block below implements the chosen method AND the alternative on identical
inputs, and reports the discrepancy. Where the alternative is not merely worse
but actively misleading, that is the point: several of these were run the wrong
way first, and the number here is what exposed it.

Decisions tested:
  D1  cross-plate vs pooled replicate pairs         (plate confound)
  D2  replicate covariance vs residual variance     (noise cancellation)
  D3  leave-one-out vs in-sample additive prior     (target leakage)
  D4  our CDI vs published cross-context r          (reproducibility confound)
  D5  Hub vs Tahoe-name mechanism annotation        (coverage)
  D6  within-compound vs pooled dose regression     (Simpson's paradox risk)

Outputs: results/tables/methodology_evidence.csv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
MIN_PAIRS = 10
N_RESP = 2000


def cor(a, b):
    a = a - a.mean(); b = b - b.mean()
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    args = ap.parse_args()
    rows = []

    from perturbmodel.evaluation.delta_eval import (build_deltas,
                                                    load_pseudobulk,
                                                    responsive_genes)
    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond, keep_plate=True)
    resp = responsive_genes(DELTA, np.ones(len(G), bool))
    D = DELTA[:, resp].astype(np.float32)
    del DELTA
    K = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                      "drug": G.drug, "conc": G.conc, "plate": G.plate})
    print(f"{len(K)} conditions, {len(resp)} response genes, "
          f"{K.line.nunique()} lines, {K.drug.nunique()} drugs", flush=True)

    # ---- D1/D2: three estimators of the same interaction, one dataset ----
    shared_t, cov_x, n_x, cov_p, n_p, var_r, n_r = 0.0, 0.0, 0, 0.0, 0, 0.0, 0
    within_r, cross_r = [], []
    for drug, gd in K.groupby("drug", observed=True):
        for _, gc in gd.groupby("conc", observed=True):
            ii = gc.i.to_numpy()
            if len(ii) < 3:
                continue
            # Leave-one-CONTEXT-out. Two weaker priors were tried first and
            # both produced a NEGATIVE pooled covariance -- the in-sample mean
            # gives about -sigma^2/n, leave-one-condition-out about
            # -2 sigma^2/(n-1), because the covariance is between two plates of
            # the same line and those priors leave the sibling plate in the
            # mean. The negative numbers were arithmetic, not absence of signal.
            ln = gc.line.to_numpy()
            if len(np.unique(ln)) < 2:
                continue
            tot = D[ii].sum(0)
            csum = {c: D[ii[ln == c]].sum(0) for c in np.unique(ln)}
            ccnt = {c: int((ln == c).sum()) for c in np.unique(ln)}
            loo, ok = {}, []
            for i, c in zip(ii, ln):
                n_out = len(ii) - ccnt[c]
                if n_out < 1:
                    continue
                loo[i] = (tot - csum[c]) / n_out; ok.append(i)
            if not ok:
                continue
            shared_t += float(np.sum([np.mean(loo[i] ** 2) for i in ok]))
            res = {i: D[i] - loo[i] for i in ok}
            ii = np.array(ok)
            for i in ii:
                var_r += float(np.mean(res[i] ** 2)); n_r += 1
            for _, gl in gc.groupby("line", observed=True):
                v = gl.i.to_numpy(); pl = gl.plate.to_numpy()
                for a in range(len(v)):
                    for b in range(a + 1, len(v)):
                        c = float(np.mean(res[v[a]] * res[v[b]]))
                        cov_p += c; n_p += 1                 # pooled (wrong)
                        r_ = cor(D[v[a]], D[v[b]])
                        if pl[a] != pl[b]:
                            cov_x += c; n_x += 1             # cross-plate only
                            cross_r.append(r_)
                        else:
                            within_r.append(r_)
    shared = shared_t / len(K)
    i_cross = cov_x / max(n_x, 1)
    i_pool = cov_p / max(n_p, 1)
    i_var = var_r / max(n_r, 1)
    sh = lambda x: max(x, 0.0) / (shared + max(x, 0.0))
    if len(within_r) == 0:
        # this build has one plate per (line, drug, dose), so no same-plate
        # replicate pair exists to compare against. The plate confound was
        # demonstrated on the dev47 build, which does carry them; reporting a
        # number here would be reporting the absence of the comparison.
        print(f"\nD1 plate handling  NOT COMPUTABLE on this build: "
              f"{len(cross_r)} cross-plate pairs and 0 within-plate pairs "
              f"(one plate per condition).\n   Cross-plate residual agreement "
              f"is r={np.median(cross_r):+.3f}. The 7x within/cross gap is "
              f"measured on the dev47 build (RESULTS.md §4).")
        rows.append({"decision": "D1 cross-plate vs pooled pairs",
                     "chosen": "cross-plate", "chosen_value": np.nan,
                     "alternative": "pool all pairs", "alt_value": np.nan,
                     "evidence": f"not computable here: 0 within-plate pairs, "
                                 f"{len(cross_r)} cross-plate"})
    else:
        ratio = np.median(within_r) / max(np.median(cross_r), 1e-9)
        print(f"\nD1 plate handling  within-plate r={np.median(within_r):+.3f} "
              f"(n={len(within_r)}) vs cross-plate r={np.median(cross_r):+.3f} "
              f"(n={len(cross_r)}), ratio {ratio:.1f}x")
        print(f"   interaction share: cross-plate only {sh(i_cross):.1%} vs "
              f"pooled {sh(i_pool):.1%}")
        rows.append({"decision": "D1 cross-plate vs pooled pairs",
                     "chosen": "cross-plate", "chosen_value": round(sh(i_cross), 4),
                     "alternative": "pool all pairs",
                     "alt_value": round(sh(i_pool), 4),
                     "evidence": f"within/cross agreement ratio {ratio:.1f}x"})
    print(f"D2 estimator       replicate covariance {i_cross:+.5f} vs residual "
          f"variance {i_var:+.5f} (raw, against shared {shared:.5f})")
    print(f"   as shares: {sh(i_cross):.1%} vs {sh(i_var):.1%} -- the variance "
          f"estimator absorbs the per-condition noise that the covariance "
          f"cancels,\n   which is the entire difference between them.")
    rows.append({"decision": "D2 covariance vs variance",
                 "chosen": "replicate covariance", "chosen_value": round(sh(i_cross), 4),
                 "alternative": "residual variance", "alt_value": round(sh(i_var), 4),
                 "evidence": f"raw {i_cross:+.5f} vs {i_var:+.5f}, "
                             f"shared {shared:.5f}"})

    # ---- D3: leave-one-out vs in-sample additive prior, at small and large n ----
    lines = sorted(K.line.unique())
    rng = np.random.default_rng(0)
    for n_ctx in (5, len(lines)):
        sub = lines if n_ctx == len(lines) else list(rng.choice(lines, 5, False))
        S = K[K.line.isin(sub)]
        loo_r, ins_r = [], []
        for (drug, conc), g in S.groupby(["drug", "conc"], observed=True):
            ii = g.i.to_numpy()
            if len(ii) < 3:
                continue
            tot = D[ii].sum(0)
            for j, i in enumerate(ii):
                loo_r.append(cor(D[i], (tot - D[i]) / (len(ii) - 1)))
                ins_r.append(cor(D[i], tot / len(ii)))
        rows.append({"decision": f"D3 additive prior at n={n_ctx} contexts",
                     "chosen": "leave-one-out", "chosen_value":
                     round(float(np.median(loo_r)), 4),
                     "alternative": "in-sample mean",
                     "alt_value": round(float(np.median(ins_r)), 4),
                     "evidence": f"{len(loo_r)} conditions"})
        print(f"D3 prior at n={n_ctx:2d}   leave-one-out r="
              f"{np.median(loo_r):+.3f}  vs in-sample r={np.median(ins_r):+.3f}"
              f"  -> leakage inflates by {np.median(ins_r)-np.median(loo_r):+.3f}")

    # ---- D4: our metric vs the published one ----
    f = TAB / "context_metric_comparison.csv"
    if f.exists():
        M = pd.read_csv(f).dropna(subset=["cdi", "r_cross", "r_replicate"])
        rc = stats.spearmanr(M.r_cross, M.r_replicate)
        print(f"D4 metric choice   published cross-context r correlates with a "
              f"compound's own\n   reproducibility at rho={rc.statistic:+.3f} "
              f"(n={len(M)}) -- it partly measures potency, not transfer")
        rows.append({"decision": "D4 CDI vs published cross-context r",
                     "chosen": "CDI (replicate-validated)", "chosen_value": np.nan,
                     "alternative": "cross-context correlation",
                     "alt_value": round(float(rc.statistic), 4),
                     "evidence": "alt confounds transfer with reproducibility"})

    # ---- D5: annotation source coverage ----
    dm = pd.read_parquet(ROOT / "data/metadata/metadata/drug_metadata.parquet")
    hub = TAB.parent.parent / "data/external/lincs/repurposing_drugs.txt"
    if hub.exists():
        import re
        nrm = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())
        H = pd.read_csv(hub, sep="\t", comment="!", low_memory=False)
        lm = pd.read_csv(ROOT / "data/external/lincs/inst_info.txt.gz", sep="\t",
                         low_memory=False)
        names = set(lm[lm.pert_type == "trt_cp"].pert_iname.map(nrm))
        c_hub = len(names & set(H.pert_iname.map(nrm)))
        c_tah = len(names & set(dm.drug.map(nrm)))
        print(f"D5 MOA source      Repurposing Hub annotates {c_hub} of "
              f"{len(names)} LINCS compounds; Tahoe names annotate {c_tah}")
        rows.append({"decision": "D5 mechanism annotation source",
                     "chosen": "Broad Repurposing Hub", "chosen_value": c_hub,
                     "alternative": "Tahoe drug names", "alt_value": c_tah,
                     "evidence": f"{len(names)} LINCS compounds"})

    # ---- D6: within-compound vs pooled dose regression ----
    f = TAB / "dose_vs_context_prism.csv"
    if f.exists():
        P = pd.read_csv(f)
        pooled = stats.spearmanr(P.dose_pos, P.cdi, nan_policy="omit")
        wr = []
        for _, g in P.groupby("compound", observed=True):
            g = g.dropna(subset=["cdi"])
            if g.dose_pos.nunique() >= 3:
                r = stats.spearmanr(g.dose_pos, g.cdi).statistic
                if np.isfinite(r):
                    wr.append(r)
        print(f"D6 dose model      within-compound median rho="
              f"{np.median(wr):+.3f} (n={len(wr)}) vs pooled across all "
              f"compounds rho={pooled.statistic:+.3f}")
        rows.append({"decision": "D6 dose trend estimation",
                     "chosen": "within-compound then sign test",
                     "chosen_value": round(float(np.median(wr)), 4),
                     "alternative": "pooled regression",
                     "alt_value": round(float(pooled.statistic), 4),
                     "evidence": f"{len(wr)} compounds"})

    E = pd.DataFrame(rows)
    E.to_csv(TAB / "methodology_evidence.csv", index=False)
    print(f"\n{len(E)} decisions with measured evidence -> "
          f"{TAB/'methodology_evidence.csv'}")
    print(E.to_string(index=False))


if __name__ == "__main__":
    main()
