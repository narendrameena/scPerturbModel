#!/usr/bin/env python3
"""Does the spectrum method beat the published index, on data with a known answer?

RESULTS.md sec.35 and sec.36 showed that an atlas-wide context-dependence index
reports near-zero when the interaction is concentrated, and that both nulls in
this project were of that kind. ``perturbmodel.spectrum`` proposes a replacement:
the published index is the TRACE of the cross-replicate covariance, so decompose
that matrix and keep the directions that survive a permuted-pairing null.

A method is only better if it is better on data whose answer is known, so this
script benchmarks it three ways before it is applied to anything.

  1. SIMULATION, varying concentration. A fixed amount of true reproducible
     interaction is placed in a subspace of r directions out of 2,000, with r
     swept from 2,000 (spread evenly) down to 2 (highly concentrated). The
     truth is identical at every point; only its arrangement changes. A good
     estimator returns the same number throughout. The trace is expected to
     hold up when the signal is spread and to collapse when it concentrates,
     because concentration is exactly when noise eigenvalues cancel it.

  2. SIMULATION, true null. No interaction at all. Any method that reports
     something here is worse than useless, and the spectrum method is the one at
     risk, because selecting the largest eigenvalues is a selection procedure --
     which is why the null edge comes from permuting the replicate pairing
     rather than from theory.

  3. THE PUBLISHED ALTERNATIVES on the same simulated data: the naive residual
     variance ratio, and ANOVA variance components (Searle et al. 1992, the
     balanced-design limit of variancePartition, Hoffman & Schadt 2016).

Only after all three does the method touch real data (``spectrum_tahoe.py``).

Outputs: results/tables/spectrum_benchmark.csv
         figure bundle results/figures/00_manuscript/spectrum_benchmark/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.spectrum import interaction_spectrum
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def simulate(n_cond=300, n_feat=2000, rank=50, var_true=0.20, sigma=1.0,
             seed=0):
    """Two replicate measurements with a known reproducible interaction.

    The interaction occupies ``rank`` orthogonal directions and carries total
    variance ``var_true`` per condition however many directions that is, so the
    truth is constant while its concentration varies. Noise is independent
    between the two replicates, which is the property every estimator here
    relies on.
    """
    rng = np.random.default_rng(seed)
    Q = np.linalg.qr(rng.normal(size=(n_feat, rank)))[0]
    scale = np.sqrt(var_true * n_feat / rank)
    S = rng.normal(0, scale, (n_cond, rank))
    G = S @ Q.T                                    # the shared interaction
    A = G + rng.normal(0, sigma, (n_cond, n_feat))
    B = G + rng.normal(0, sigma, (n_cond, n_feat))
    truth = float((G * G).sum() / (n_cond * n_feat))
    return A, B, truth


def trace_index(A, B):
    """The published-style index: mean cross-replicate inner product."""
    return float(np.mean(np.sum(A * B, axis=1)) / A.shape[1])


def variance_index(A, B):
    """Residual variance rather than covariance -- the naive alternative."""
    return float(np.mean((0.5 * (A + B)) ** 2))


def anova_index(A, B):
    """ANOVA/Henderson components on the same data, as two 'replicates'.

    Balanced two-level design: the condition is the group, the two replicates
    are the observations. sigma_between estimates the reproducible variance,
    which is the interaction here since no other structure was simulated.
    """
    Y = np.stack([A, B], axis=1)                  # cond x rep x feat
    cell = Y.mean(axis=1)
    gm = cell.mean(axis=0)
    n_rep = 2
    ms_b = (n_rep * ((cell - gm) ** 2).sum(0) / max(len(cell) - 1, 1))
    ms_w = ((Y - cell[:, None, :]) ** 2).sum((0, 1)) / max(len(cell), 1)
    v_b = np.maximum((ms_b - ms_w) / n_rep, 0)
    return float(v_b.mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=120)
    ap.add_argument("--n-seeds", type=int, default=3)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    ranks = [2000, 800, 300, 100, 30, 10, 4, 2]
    rows = []
    print("1. CONCENTRATION SWEEP — the truth is constant at every rank",
          flush=True)
    for r in ranks:
        for seed in range(args.n_seeds):
            A, B, truth = simulate(rank=r, seed=seed)
            sp = interaction_spectrum(A, B, n_perm=args.n_perm, seed=seed)
            rows.append({"setting": "concentration", "rank": r, "seed": seed,
                         "truth": truth, "trace": sp.trace,
                         "spectrum": sp.reproducible,
                         "n_components": sp.n_components,
                         "variance": variance_index(A, B),
                         "anova": anova_index(A, B)})
        g = pd.DataFrame(rows)
        g = g[(g.setting == "concentration") & (g["rank"] == r)]
        print(f"   rank {r:5d}: truth {g.truth.mean():.4f}   "
              f"trace {g.trace.mean():+.4f}   spectrum "
              f"{g.spectrum.mean():+.4f}   ({g.n_components.mean():.0f} "
              f"components)", flush=True)

    print("\n2. TRUE NULL — no interaction at all", flush=True)
    for seed in range(args.n_seeds):
        rng = np.random.default_rng(100 + seed)
        A = rng.normal(0, 1.0, (300, 2000))
        B = rng.normal(0, 1.0, (300, 2000))
        sp = interaction_spectrum(A, B, n_perm=args.n_perm, seed=seed)
        rows.append({"setting": "null", "rank": 0, "seed": seed, "truth": 0.0,
                     "trace": sp.trace, "spectrum": sp.reproducible,
                     "n_components": sp.n_components,
                     "variance": variance_index(A, B),
                     "anova": anova_index(A, B)})
    N = pd.DataFrame(rows)
    N = N[N.setting == "null"]
    print(f"   trace {N.trace.mean():+.5f}   spectrum "
          f"{N.spectrum.mean():+.5f}   components "
          f"{N.n_components.mean():.1f}   variance index "
          f"{N.variance.mean():.3f}   ANOVA {N.anova.mean():+.5f}")
    print("   The spectrum method must report ~0 here; selecting the largest "
          "eigenvalues\n   is a selection procedure, and the permuted-pairing "
          "edge is what controls it.")

    T = pd.DataFrame(rows)
    T.to_csv(TAB / "spectrum_benchmark.csv", index=False)

    C = T[T.setting == "concentration"].groupby("rank").mean(numeric_only=True)
    C = C.sort_index(ascending=False)
    print("\n3. ACCURACY over the sweep (relative error against the truth)")
    for col in ("trace", "spectrum", "anova"):
        err = float(np.mean(np.abs(C[col] - C.truth) / C.truth))
        worst = float(np.max(np.abs(C[col] - C.truth) / C.truth))
        print(f"   {col:9s} mean |error| {err:6.1%}   worst {worst:6.1%}")
    lo = C.index.min()
    print(f"\n   at the most concentrated setting (rank {lo}): "
          f"truth {C.truth.loc[lo]:.4f}, trace {C.trace.loc[lo]:+.4f} "
          f"({C.trace.loc[lo]/C.truth.loc[lo]:.0%} of it), "
          f"spectrum {C.spectrum.loc[lo]:+.4f} "
          f"({C.spectrum.loc[lo]/C.truth.loc[lo]:.0%})")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)
    x = C.index.to_numpy()
    ax[0].plot(x, C.truth, "k--", lw=1.6, label="truth (constant)")
    ax[0].plot(x, C.trace, "o-", color=ORANGE, lw=2, ms=6,
               label="trace (published-style index)")
    ax[0].plot(x, C.spectrum, "o-", color=VIOLET, lw=2, ms=6,
               label="spectrum (this work)")
    ax[0].plot(x, C.anova, "o-", color=AQUA, lw=1.6, ms=5,
               label="ANOVA components")
    ax[0].set_xscale("log")
    ax[0].invert_xaxis()
    ax[0].set_xlabel("directions the interaction occupies (of 2,000)")
    ax[0].set_ylabel("estimated interaction variance")
    ax[0].legend(frameon=False, fontsize=7)
    ax[0].set_title("a  Same truth, varying concentration", loc="left",
                    fontweight="bold", fontsize=9.5)

    rel = (C[["trace", "spectrum", "anova"]].div(C.truth, axis=0))
    ax[1].axhline(1.0, ls="--", color="#555", lw=1.4)
    for col, c in (("trace", ORANGE), ("spectrum", VIOLET), ("anova", AQUA)):
        ax[1].plot(x, rel[col], "o-", color=c, lw=2, ms=6, label=col)
    ax[1].set_xscale("log"); ax[1].invert_xaxis()
    ax[1].set_xlabel("directions the interaction occupies")
    ax[1].set_ylabel("estimate / truth")
    ax[1].legend(frameon=False, fontsize=7)
    ax[1].set_title("b  Recovery, where 1.0 is correct", loc="left",
                    fontweight="bold", fontsize=9.5)

    A, B, truth = simulate(rank=10, seed=0)
    sp = interaction_spectrum(A, B, n_perm=args.n_perm)
    k = min(60, len(sp.eigvals))
    ax[2].plot(np.arange(1, k + 1), sp.eigvals[:k], "o-", color=VIOLET, lw=1.6,
               ms=4, label="observed")
    ax[2].axhline(sp.null_edge, ls="--", color=ORANGE, lw=1.6,
                  label="permuted-pairing edge")
    ax[2].axhline(0, color="#888", lw=0.9)
    ax[2].set_xlabel("component")
    ax[2].set_ylabel("eigenvalue")
    ax[2].legend(frameon=False, fontsize=7.5)
    ax[2].text(0.5, 0.55, f"trace = {sp.trace:+.4f}\nsum of the "
               f"{sp.n_components} components\nabove the edge = "
               f"{sp.reproducible:+.4f}", transform=ax[2].transAxes,
               ha="center", fontsize=7.5, color="#444")
    ax[2].set_title("c  Why the trace loses it (rank 10)", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("The published index is a trace; the interaction it misses is "
                 "in the spectrum", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    d = save_figure(fig, "spectrum_benchmark", FIG,
                    source_data={"sweep": T}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
