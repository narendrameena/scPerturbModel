#!/usr/bin/env python3
"""Which drug mechanisms are most context-dependent, and is that predictable?

Everything so far has asked what makes a CELL LINE respond differently. This
asks the complementary, positive-framed question: which DRUGS have effects that
are conserved across cellular contexts, and which have effects that are mostly
context-specific — and can that be predicted from chemistry or target class
before running the experiment?

Per drug d, over its (line, dose, plate) conditions and the responsive genes:

  conserved(d)   mean square of the drug's own dose-matched mean response
                 (what the drug does in every line)
  context(d)     REPRODUCIBLE line-specific variance, estimated as the
                 covariance of residuals between different doses of the same
                 (line, drug) measured on DIFFERENT plates — independent noise
                 and plate batch structure both cancel
  CDI(d)         context / (conserved + context)
                 0 = the drug does the same thing everywhere
                 1 = the drug's reproducible effect is entirely line-specific

Then: rank drugs, compare mechanism classes (Kruskal-Wallis), and test whether
CDI is predictable from chemical structure alone via cross-validated ridge
regression on ECFP4 fingerprints — a prospective question, since fingerprints
are available before any experiment.

Outputs: results/tables/drug_context_dependence.csv
         figure bundle results/figures/10_drugs/drug_context_dependence/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.evaluation.delta_eval import (build_deltas, load_pseudobulk,
                                                responsive_genes)
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures" / "10_drugs"
TAB = ROOT / "results" / "tables"
MIN_LINES = 15
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    ap.add_argument("--tag", default="full")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond, keep_plate=True)
    resp = responsive_genes(DELTA, np.ones(len(G), bool))
    D = DELTA[:, resp]
    print(f"{len(G)} plate-level conditions, {len(resp)} responsive genes",
          flush=True)

    meta = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                         "drug": G.drug, "conc": G.conc, "plate": G.plate})
    rng = np.random.default_rng(0)
    recs = []
    for drug, gd in meta.groupby("drug", observed=True):
        if gd.line.nunique() < MIN_LINES:
            continue
        # drug's own conserved effect: leave-one-LINE-out mean per dose. An
        # in-sample mean was used here originally and is wrong: it makes the
        # residuals of a dose group sum to zero, biasing every residual pair
        # negative.
        conserved = 0.0
        resid = np.zeros((len(gd), D.shape[1]), dtype=np.float32)
        pos = {v: k for k, v in enumerate(gd.i.to_numpy())}
        have = set()
        for cc, gc in gd.groupby("conc", observed=True):
            idx = gc.i.to_numpy(); ln_ = gc.line.to_numpy()
            if len(np.unique(ln_)) < 2:
                continue
            tot = D[idx].sum(0)
            csum = {c: D[idx[ln_ == c]].sum(0) for c in np.unique(ln_)}
            ccnt = {c: int((ln_ == c).sum()) for c in np.unique(ln_)}
            for j, c in zip(idx, ln_):
                n_out = len(idx) - ccnt[c]
                if n_out < 1:
                    continue
                loo = (tot - csum[c]) / n_out
                conserved += float(np.mean(loo ** 2))
                resid[pos[j]] = D[j] - loo
                have.add(j)
        if not have:
            continue
        conserved /= len(have)

        # reproducible line-specific variance: cross-plate pairs of the SAME
        # line, minus a matched cross-LINE null. Both carry the residual-
        # construction offset; only the first carries interaction, so the
        # difference is the interaction. Without this subtraction the estimate
        # is the offset plus the signal and comes out negative for most drugs.
        cov, npair = 0.0, 0
        for (ln,), gl in gd.groupby(["line"], observed=True):
            v = [x for x in gl.i.to_numpy() if x in have]
            pl = {x: p_ for x, p_ in zip(gl.i.to_numpy(), gl.plate.to_numpy())}
            for a in range(len(v)):
                for b in range(a + 1, len(v)):
                    if pl[v[a]] == pl[v[b]]:
                        continue
                    cov += float(np.mean(resid[pos[v[a]]] * resid[pos[v[b]]]))
                    npair += 1
        if npair < 10:
            continue
        keys = np.array(sorted(have))
        kl = gd.set_index("i").line.loc[keys].to_numpy()
        kp = gd.set_index("i").plate.loc[keys].to_numpy()
        c3, n3 = 0.0, 0
        for _ in range(min(len(keys) * 6, 3000)):
            a, b = rng.integers(0, len(keys), 2)
            if kl[a] == kl[b] or kp[a] == kp[b]:
                continue
            c3 += float(np.mean(resid[pos[keys[a]]] * resid[pos[keys[b]]]))
            n3 += 1
        off = (c3 / n3) if n3 else 0.0
        context = cov / npair - off
        denom = conserved + max(context, 0.0)
        recs.append({"drug": drug, "n_lines": gd.line.nunique(),
                     "n_conditions": len(gd), "n_pairs": npair,
                     "n_null_pairs": n3, "raw_cov": cov / npair,
                     "null_offset": off,
                     "conserved": conserved, "context": context,
                     "cdi": float(max(context, 0.0) / denom) if denom > 0 else np.nan,
                     "total_effect": conserved + max(context, 0.0)})
    res = pd.DataFrame(recs).dropna(subset=["cdi"])
    dm = pd.read_parquet(ROOT / "data/metadata/metadata/drug_metadata.parquet")
    res = res.merge(dm[["drug", "moa-fine", "targets", "human-approved"]],
                    on="drug", how="left")
    res = res.sort_values("cdi", ascending=False)
    res.to_csv(TAB / "drug_context_dependence.csv", index=False)
    print(f"\n{len(res)} drugs scored. CDI median {res.cdi.median():.3f}")
    print("\nMOST context-dependent drugs:")
    print(res.head(12)[["drug", "moa-fine", "cdi", "total_effect"]]
          .to_string(index=False))
    print("\nMOST conserved drugs (lowest CDI):")
    print(res.tail(10)[["drug", "moa-fine", "cdi", "total_effect"]]
          .to_string(index=False))

    # ---------------- by mechanism ----------------
    known = res[res["moa-fine"].notna() & (res["moa-fine"] != "unclear")]
    grp = (known.groupby("moa-fine").cdi
           .agg(["median", "size"]).query("size >= 3")
           .sort_values("median", ascending=False))
    print(f"\ncontext-dependence by mechanism (n>=3 drugs):")
    print(grp.round(3).to_string())
    if len(grp) > 2:
        samples = [known[known["moa-fine"] == m].cdi.to_numpy()
                   for m in grp.index]
        h, p = stats.kruskal(*samples)
        print(f"Kruskal-Wallis across {len(grp)} mechanisms: H={h:.1f}, p={p:.2e}")

    # ---------------- predictable from chemistry? ----------------
    ecfp = np.load(ROOT / "data/processed/drug_ecfp.npz", allow_pickle=True)
    fp_of = dict(zip(ecfp["drugs"].tolist(), ecfp["fp"]))
    sub = res[res.drug.isin(fp_of)]
    Fp = np.stack([fp_of[d] for d in sub.drug]).astype(np.float32)
    y = sub.cdi.to_numpy()
    rng = np.random.default_rng(0)
    folds = rng.permutation(len(y)) % 5
    pred = np.zeros_like(y)
    for k in range(5):
        tr, te = folds != k, folds == k
        A = Fp[tr]
        A = A - A.mean(0)
        yt = y[tr] - y[tr].mean()
        w = np.linalg.solve(A.T @ A + 50.0 * np.eye(A.shape[1]), A.T @ yt)
        pred[te] = (Fp[te] - Fp[tr].mean(0)) @ w + y[tr].mean()
    r_chem = float(np.corrcoef(pred, y)[0, 1])
    print(f"\nchemistry -> context-dependence: cross-validated r = {r_chem:.3f} "
          f"(n={len(y)} drugs, 5-fold ridge on ECFP4)")

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2), constrained_layout=True)

    top = grp.head(14)
    yy = np.arange(len(top))[::-1]
    axes[0].barh(yy, top["median"], color=BLUE, height=0.65)
    axes[0].set_yticks(yy, [f"{m[:34]} (n={int(n)})"
                            for m, n in zip(top.index, top["size"])], fontsize=7.5)
    axes[0].axvline(res.cdi.median(), color="#888888", ls="--", lw=0.9)
    axes[0].set_xlabel("context-dependence index (median)")
    axes[0].set_title("A  Which mechanisms are context-dependent?",
                      loc="left", fontweight="bold", fontsize=10)

    axes[1].scatter(res.total_effect, res.cdi, s=10, alpha=0.45, color=VIOLET,
                    edgecolors="none")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("total reproducible effect size")
    axes[1].set_ylabel("context-dependence index")
    rr = stats.spearmanr(res.total_effect, res.cdi)
    axes[1].set_title(f"B  Potency vs context (rho={rr.statistic:+.2f})",
                      loc="left", fontweight="bold", fontsize=10)

    axes[2].scatter(pred, y, s=10, alpha=0.45, color=AQUA, edgecolors="none")
    lim = [min(pred.min(), y.min()), max(pred.max(), y.max())]
    axes[2].plot(lim, lim, ls="--", color="#888888", lw=1)
    axes[2].set_xlabel("predicted CDI from ECFP4 (5-fold CV)")
    axes[2].set_ylabel("observed CDI")
    axes[2].set_title(f"C  Predictable from chemistry? r={r_chem:.2f}",
                      loc="left", fontweight="bold", fontsize=10)

    fig.suptitle("Context-dependence is a property of the drug: which "
                 "mechanisms transfer across cell lines?",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "drug_context_dependence", FIG,
                    source_data={"drugs": res,
                                 "by_mechanism": grp.reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
