#!/usr/bin/env python3
"""Does sparse detection beat the pooled test where the pooled test fails?

RESULTS.md sec.38 rejected a proposed replacement ESTIMATOR: the trace recovers
planted interaction to within 0.5% at every concentration, so there is nothing
wrong with it as a measure of magnitude. What failed three times in this project
was the TEST -- asking whether a pooled variance fraction differs from zero, when
the alternative is a handful of strongly interacting pairs among thousands of
null ones. That is the classical sparse-detection regime, where a sum has
vanishing power (Donoho & Jin, Ann. Statist. 2004).

``perturbmodel.sparse_interaction`` replaces the test, not the estimator:
per-pair cross-replicate reproducibility, FDR over pairs, Higher Criticism for
global detection, and a greedy search for the reproducible submatrix.

The benchmark is a power curve. A fixed total amount of interaction is spread
over ``k`` of 400 conditions, with ``k`` swept from 400 (dense) to 4 (sparse).
Total signal is held constant, so any change in power is a property of the test
rather than of the data's information content. Three tests see identical data:

  pooled          permutation test on the trace -- what the field does
  FDR over pairs  count of pairs reproducing at q < 0.05
  Higher Criticism  global sparse-detection statistic

and a true null is run to confirm none of them fires when nothing is there --
the risk being that per-pair testing and a maximum statistic are both selection
procedures.

Finally the biclustering is tested on a planted block, since the biology that
motivated it (one drug class acting on one genotype group) is a submatrix rather
than scattered cells.

Outputs: results/tables/sparse_benchmark.csv
         figure bundle results/figures/00_manuscript/sparse_benchmark/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.sparse_interaction import detect, find_bicluster
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def simulate_sparse(n_cond=400, n_feat=300, k=40, total=300.0, sigma=1.0,
                    seed=0):
    """Interaction in k of n conditions, with the TOTAL held constant.

    Each affected condition gets an effect of size sqrt(total / k), so the sum
    of squared effects is ``total`` however sparse the arrangement. A test whose
    power depends only on total signal would be flat across k; a sum-based test
    is not.

    ``total`` is set so both regimes are visible. At k = 400 each condition
    carries a cosine of about 0.003 against a null spread of 1/sqrt(300) = 0.058,
    invisible per pair but easily summed; at k = 4 each carries about 0.25, a
    four-sigma per-pair effect that the sum cannot see. A first version used
    total = 40, which put every point below detection for every test and
    produced an all-zero table.
    """
    rng = np.random.default_rng(seed)
    A = rng.normal(0, sigma, (n_cond, n_feat))
    B = rng.normal(0, sigma, (n_cond, n_feat))
    idx = rng.choice(n_cond, k, replace=False) if k else np.array([], int)
    amp = np.sqrt(total / max(k, 1))
    for i in idx:
        v = rng.normal(size=n_feat)
        v /= np.linalg.norm(v)
        A[i] += amp * v
        B[i] += amp * v
    return A, B, set(idx.tolist())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=300)
    ap.add_argument("--n-seeds", type=int, default=8)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    ks = [400, 200, 100, 40, 16, 8, 4]
    rows = []
    print("1. POWER CURVE — total interaction held constant, spread over k of "
          "400 conditions", flush=True)
    for k in ks:
        for seed in range(args.n_seeds):
            A, B, truth = simulate_sparse(k=k, seed=seed)
            r = detect(A, B, n_perm=args.n_perm, seed=seed)
            hit = (r.table.q < 0.05).to_numpy()
            tp = len(set(np.where(hit)[0]) & truth)
            rows.append({"setting": "power", "k": k, "seed": seed,
                         "pooled_p": r.trace_p, "hc_p": r.hc_p,
                         "n_sig": r.n_sig, "true_positives": tp,
                         "recall": tp / max(len(truth), 1),
                         "precision": tp / max(r.n_sig, 1)})
        g = pd.DataFrame(rows)
        g = g[(g.setting == "power") & (g.k == k)]
        print(f"   k={k:4d}: pooled rejects {(g.pooled_p<0.05).mean():5.0%}   "
              f"HC rejects {(g.hc_p<0.05).mean():5.0%}   "
              f"FDR finds {g.n_sig.mean():5.1f} pairs "
              f"(recall {g.recall.mean():.0%}, "
              f"precision {g.precision.mean():.0%})", flush=True)

    print("\n2. TRUE NULL — no interaction", flush=True)
    for seed in range(args.n_seeds):
        A, B, _ = simulate_sparse(k=0, total=0.0, seed=200 + seed)
        r = detect(A, B, n_perm=args.n_perm, seed=seed)
        rows.append({"setting": "null", "k": 0, "seed": seed,
                     "pooled_p": r.trace_p, "hc_p": r.hc_p, "n_sig": r.n_sig,
                     "true_positives": 0, "recall": np.nan,
                     "precision": np.nan})
    N = pd.DataFrame(rows)
    N = N[N.setting == "null"]
    print(f"   pooled false-positive rate {(N.pooled_p<0.05).mean():.0%}   "
          f"HC {(N.hc_p<0.05).mean():.0%}   "
          f"FDR pairs called {N.n_sig.mean():.1f}")
    print("   Both new tests are selection procedures, so a clean null here is "
          "what makes\n   the power curve above meaningful.")

    T = pd.DataFrame(rows)
    T.to_csv(TAB / "sparse_benchmark.csv", index=False)

    print("\n3. BICLUSTER — a drug class acting on a genotype group is a block",
          flush=True)
    rng = np.random.default_rng(0)
    n_ctx, n_pert = 40, 60
    Z = rng.normal(0, 1, (n_ctx, n_pert))
    br, bc = sorted(rng.choice(n_ctx, 8, replace=False)), \
        sorted(rng.choice(n_pert, 5, replace=False))
    Z[np.ix_(br, bc)] += 1.6
    res = find_bicluster(Z, n_perm=400)
    rec = len(set(res["rows"]) & set(br)) / len(br)
    pre = len(set(res["rows"]) & set(br)) / max(len(res["rows"]), 1)
    recc = len(set(res["cols"]) & set(bc)) / len(bc)
    print(f"   planted {len(br)} contexts x {len(bc)} perturbations")
    print(f"   found   {res['n_rows']} x {res['n_cols']}   "
          f"score {res['score']:.2f} vs null {res['null_mean']:.2f}, "
          f"p = {res['p']:.4f}")
    print(f"   context recall {rec:.0%}, precision {pre:.0%}; "
          f"perturbation recall {recc:.0%}")

    P = T[T.setting == "power"].groupby("k").mean(numeric_only=True)
    P = P.sort_index(ascending=False)
    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)
    pw = T[T.setting == "power"].groupby("k").agg(
        pooled=("pooled_p", lambda s: float((s < 0.05).mean())),
        hc=("hc_p", lambda s: float((s < 0.05).mean()))).sort_index(
            ascending=False)
    x = pw.index.to_numpy()
    ax[0].plot(x, pw.pooled, "o-", color=ORANGE, lw=2, ms=7,
               label="pooled test (the field's)")
    ax[0].plot(x, pw.hc, "o-", color=VIOLET, lw=2, ms=7,
               label="Higher Criticism (sparse)")
    ax[0].axhline(0.05, ls=":", color="#888", lw=1.2)
    ax[0].set_xscale("log"); ax[0].invert_xaxis()
    ax[0].set_xlabel("conditions carrying the interaction (of 400)")
    ax[0].set_ylabel("power at α = 0.05")
    ax[0].set_ylim(-0.03, 1.05)
    ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_title("a  Same total signal, varying sparsity", loc="left",
                    fontweight="bold", fontsize=9.5)

    ax[1].plot(x, P.recall, "o-", color=AQUA, lw=2, ms=7, label="recall")
    ax[1].plot(x, P.precision, "o-", color=BLUE, lw=2, ms=7, label="precision")
    ax[1].set_xscale("log"); ax[1].invert_xaxis()
    ax[1].set_xlabel("conditions carrying the interaction")
    ax[1].set_ylabel("per-pair FDR performance")
    ax[1].set_ylim(-0.03, 1.05)
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].set_title("b  And it says WHICH pairs", loc="left",
                    fontweight="bold", fontsize=9.5)

    im = ax[2].imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
    for r_ in res["rows"]:
        for c_ in res["cols"]:
            ax[2].add_patch(plt.Rectangle((c_ - .5, r_ - .5), 1, 1, fill=False,
                                          edgecolor="black", lw=0.8))
    ax[2].set_xlabel("perturbation"); ax[2].set_ylabel("context")
    ax[2].text(0.5, -0.22, f"recovered {rec:.0%} of the planted block "
               f"(p = {res['p']:.3f})", transform=ax[2].transAxes,
               ha="center", fontsize=7.5, color="#444")
    fig.colorbar(im, ax=ax[2], fraction=0.046)
    ax[2].set_title("c  Blocks, because that is the biology", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Replacing the test, not the estimator: sparse detection "
                 "where a pooled index has no power", fontsize=10.5, x=0.005,
                 ha="left", fontweight="bold")
    d = save_figure(fig, "sparse_benchmark", FIG,
                    source_data={"benchmark": T,
                                 "bicluster": pd.DataFrame([res])},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
