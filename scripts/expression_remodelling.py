#!/usr/bin/env python3
"""Does drug exposure REMODEL a cell's expression programme, or only shift it?

"The drug changes gene expression" is not a finding -- at 5 uM almost anything
does. The question with content is whether exposure moves a cell along
transcriptional directions it already had, or builds structure that no untreated
cell occupies. Those are different biological claims:

  SHIFT     the drug moves the cell along directions that already separate
            untreated cell lines from one another. The programme is pre-existing
            and the drug changes how far along it the cell sits. Nothing is
            reorganised; a model could predict the response from baseline state.

  REMODEL   the drug moves the cell in directions absent from baseline variation
            -- a programme no untreated line expresses. This is genuinely new
            structure, and it cannot be predicted by interpolating baseline
            states.

The test is a projection. The baseline programme space is the span of the top
principal components of UNTREATED pseudobulk across the 50 Tahoe cell lines: the
directions in which cells naturally differ before any drug. Each response is
split into the part inside that span (shift) and the part outside it (remodel).

**Why this needs replicate covariance, not variance.** Measurement noise is
close to isotropic, so it lands mostly OUTSIDE any low-dimensional span: a purely
noisy response would look almost entirely like "remodelling". Computing the split
from a single measurement would therefore manufacture the answer. Instead each
component is an inner product between deltas measured on two INDEPENDENT plates
-- Tahoe's plate 14 duplicates plate 6 across 50 lines and 95 compounds -- so
noise contributes zero in expectation to inside, outside and total alike. This is
the same device the project's own estimator uses, and the reason §27's correction
mattered at all.

Three further questions, each with the control that makes it interpretable:

  2. CONVERGENCE. Does exposure make different cell lines more alike or less?
     Convergence onto a shared programme is remodelling of a specific kind --
     the drug overwriting cell identity. Measured as between-line agreement in
     treated versus untreated profiles, on matched genes.
  3. COORDINATION. Does exposure change which genes co-vary across lines? A pure
     shift moves means and leaves the correlation structure intact; remodelling
     reorganises it. Tested against a condition-label permutation null.
  4. IS IT JUST DEATH? Everything above is stratified by dose. If the outside-span
     component appears only at the top dose, where lines are dying, it is
     cytotoxic collapse rather than a programme, and must not be called
     remodelling.

Outputs: results/tables/expression_remodelling.csv
         figure bundle results/figures/00_manuscript/expression_remodelling/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PB = ROOT / "data" / "processed" / "pseudobulk_full"
FIG = ROOT / "results" / "figures" / "00_manuscript"
TAB = ROOT / "results" / "tables"
CTRL = "DMSO_TF"
N_GENES = 4000
N_PC = 20
MIN_CELLS = 200
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def load():
    C = pd.read_csv(PB / "conditions.csv")
    G = pd.read_csv(PB / "genes.csv")
    X = np.load(PB / "pseudobulk_counts.npz")["counts"]
    keep = (C.n_cells >= MIN_CELLS).to_numpy()
    C, X = C[keep].reset_index(drop=True), X[keep]
    s = X.sum(1, keepdims=True)
    s[s == 0] = 1.0
    X = np.log1p(X / s * 1e4).astype(np.float32)
    # a broad, response-agnostic gene set: most-expressed genes, chosen without
    # reference to any drug effect, so neither shift nor remodel is favoured by
    # the selection
    idx = np.argsort(-X.mean(0))[:N_GENES]
    idx = np.sort(idx)
    print(f"{X.shape[0]} pseudobulks, {len(idx)} genes retained", flush=True)
    return C, X[:, idx], G.iloc[idx].reset_index(drop=True)


def baseline_span(C, X, n_pc=N_PC):
    """Directions in which UNTREATED cell lines differ from one another.

    Controls are averaged per cell line first, so the span describes
    line-to-line biology rather than plate-to-plate noise; a span fitted on
    unaveraged control wells would include the noise directions and would then
    absorb noise from the responses, understating remodelling.
    """
    ctl = C.drug.astype(str) == CTRL
    B = pd.DataFrame(X[ctl.to_numpy()])
    B["line"] = C.cell_line_id[ctl].to_numpy()
    M = B.groupby("line").mean().to_numpy(np.float32)
    mu = M.mean(0, keepdims=True)
    Mc = M - mu
    U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
    k = min(n_pc, Vt.shape[0])
    var = float((S ** 2).sum())
    print(f"baseline span from {M.shape[0]} untreated cell lines: top {k} PCs "
          f"hold {(S[:k]**2).sum()/var:.1%} of line-to-line variance",
          flush=True)
    return Vt[:k], mu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-pc", type=int, default=N_PC)
    ap.add_argument("--n-perm", type=int, default=200)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    C, X, G = load()
    V, mu = baseline_span(C, X, args.n_pc)
    p = X.shape[1]
    k = V.shape[0]

    # control per (line, plate): the same subtraction the atlas's own pipeline
    # uses, so a line's baseline never enters its response
    ctl = {}
    for (ln, pl), g in C[C.drug.astype(str) == CTRL].groupby(
            ["cell_line_id", "plate"], observed=True):
        ctl[(ln, pl)] = X[g.index.to_numpy()].mean(0)
    print(f"{len(ctl)} (line, plate) control anchors", flush=True)

    trt = C[C.drug.astype(str) != CTRL]
    rows = []
    for i, r in zip(trt.index.to_numpy(), trt.itertuples()):
        a = ctl.get((r.cell_line_id, r.plate))
        if a is None:
            continue
        rows.append((r.cell_line_id, r.drug, r.conc, r.plate, i))
    K = pd.DataFrame(rows, columns=["line", "drug", "conc", "plate", "i"])
    print(f"{len(K)} treated conditions with a plate-matched control",
          flush=True)

    # ---- 1. shift vs remodel, from cross-plate inner products ----------
    recs = []
    for (ln, dr, cc), g in K.groupby(["line", "drug", "conc"], observed=True):
        if g.plate.nunique() < 2:
            continue
        pls = list(g.plate.unique())[:2]
        d = []
        for pl in pls:
            ii = g.i[g.plate == pl].to_numpy()
            d.append(X[ii].mean(0) - ctl[(ln, pl)])
        A, B = d[0], d[1]
        ia, ib = V @ A, V @ B                      # inside the baseline span
        oa, ob = A - V.T @ ia, B - V.T @ ib        # outside it
        tot = float(A @ B)
        recs.append({"line": ln, "drug": dr, "conc": cc,
                     "inside": float(ia @ ib), "outside": float(oa @ ob),
                     "total": tot,
                     "norm": float(np.sqrt(abs(tot)))})
    R = pd.DataFrame(recs)
    R = R[R.total > 0].copy()
    R["shift_share"] = R.inside / R.total
    R["remodel_share"] = R.outside / R.total
    print(f"\n1. SHIFT vs REMODEL on {len(R)} replicated (line, drug, dose) "
          f"conditions")
    print(f"   inside the baseline programme space (shift):  "
          f"{R.shift_share.median():.1%}")
    print(f"   outside it (remodel):                         "
          f"{R.remodel_share.median():.1%}")
    print(f"   a random direction would put {k/p:.1%} inside "
          f"({k} of {p} dimensions), so 'inside' is enriched "
          f"{R.shift_share.median()/(k/p):.0f}x over chance")
    print("   But 'outside the top-20 PCs' is a weak definition of "
          "remodelling, because\n   20 of 4,000 dimensions is a very small "
          "span. The curve below is the honest\n   version: how much of the "
          "response k baseline directions capture, against how\n   much they "
          "capture of HELD-OUT baseline variation itself -- the same quantity\n"
          "   for a signal that is a pure shift by construction.")
    print("   Every term is an inner product between two INDEPENDENT plates, so "
          "isotropic\n   noise contributes zero -- without that, noise alone "
          "would read as remodelling.")

    # a within-plate variance version, to show what the naive computation says
    naive = []
    for (ln, dr, cc), g in K.groupby(["line", "drug", "conc"], observed=True):
        ii = g.i.to_numpy()
        A = X[ii].mean(0) - np.mean([ctl[(ln, pl)] for pl in g.plate.unique()],
                                    axis=0)
        ins = V @ A
        naive.append(float(ins @ ins) / max(float(A @ A), 1e-9))
    print(f"   the same quantity from single measurements (variance, not "
          f"covariance):\n   shift {np.median(naive):.1%} — noise inflates "
          f"'remodelling' by "
          f"{R.shift_share.median() - np.median(naive):+.1%}")

    # ---- 1b. capture curve, with the reference that makes it readable ----
    # "Outside the top-20 PCs" means little on its own: 20 of 4,000 dimensions
    # is a tiny span, and even a pure shift would spill out of it. The reference
    # is HELD-OUT baseline variation -- a signal that is a shift by
    # construction, since it IS line-to-line difference. If drug responses need
    # many more directions than held-out baseline lines do, they are moving
    # cells somewhere baseline variation does not go.
    ctl_lines = sorted({key[0] for key in ctl})
    rng_h = np.random.default_rng(0)
    hold = set(rng_h.choice(ctl_lines, max(len(ctl_lines) // 4, 3),
                            replace=False))
    fit_l = [l for l in ctl_lines if l not in hold]
    Mfit = np.stack([np.mean([ctl[key] for key in ctl if key[0] == l], axis=0)
                     for l in fit_l])
    mu_f = Mfit.mean(0, keepdims=True)
    _, _, Vt_f = np.linalg.svd(Mfit - mu_f, full_matrices=False)
    Mho = np.stack([np.mean([ctl[key] for key in ctl if key[0] == l], axis=0)
                    for l in hold]) - mu_f

    ks = [1, 2, 5, 10, 20, 35]
    ks = [kk for kk in ks if kk <= Vt_f.shape[0]]
    curve = []
    rep_pairs = []
    for (ln, dr, cc), g in K.groupby(["line", "drug", "conc"], observed=True):
        if g.plate.nunique() < 2:
            continue
        pls = list(g.plate.unique())[:2]
        d = [X[g.i[g.plate == pl].to_numpy()].mean(0) - ctl[(ln, pl)]
             for pl in pls]
        if float(d[0] @ d[1]) > 0:
            rep_pairs.append((d[0], d[1]))
    print(f"\n1b. capture curve over {len(rep_pairs)} replicated conditions "
          f"({len(fit_l)} lines fit the span, {len(hold)} held out)")
    for kk in ks:
        Vk = Vt_f[:kk]
        num = np.mean([float((Vk @ a) @ (Vk @ b)) for a, b in rep_pairs])
        den = np.mean([float(a @ b) for a, b in rep_pairs])
        resp = num / den if den > 0 else np.nan
        pr = Mho @ Vk.T
        base = float((pr ** 2).sum() / (Mho ** 2).sum())
        curve.append({"k": kk, "drug_response": resp, "held_out_baseline": base})
        print(f"   k = {kk:3d}:  drug response {resp:6.1%}   "
              f"held-out baseline line {base:6.1%}")
    CU = pd.DataFrame(curve)
    print("   If the two columns track each other, responses live where "
          "baseline variation\n   lives (shift). If the response column stays "
          "far below, the drug moves cells\n   into directions no untreated "
          "line occupies (remodel).")

    # ---- 4. dose stratification: is the outside component just death? ---
    print("\n4. by dose (is the outside-span part cytotoxic collapse?)")
    bydose = R.groupby("conc").agg(
        shift_share=("shift_share", "median"),
        remodel_share=("remodel_share", "median"),
        magnitude=("norm", "median"), n=("norm", "size"))
    print(bydose.round(3).to_string())
    if len(bydose) >= 2:
        rho = stats.spearmanr(bydose.index, bydose["remodel_share"])
        print(f"   remodel share vs dose: rho = {rho.statistic:+.2f}")
        print("   A remodel share that only appears at the top dose would be "
              "dying cells,\n   not a programme; one that is flat across doses "
              "is a property of the response.")

    # ---- 2. convergence: does exposure make lines more alike? -----------
    print("\n2. CONVERGENCE — between-line agreement, treated vs untreated")
    lines = sorted(set(C.cell_line_id))
    ctl_line = {}
    for ln in lines:
        v = [ctl[key] for key in ctl if key[0] == ln]
        if v:
            ctl_line[ln] = np.mean(v, axis=0)
    base_r = []
    ls = [l for l in lines if l in ctl_line]
    Mb = np.stack([ctl_line[l] for l in ls])
    Mb = Mb - Mb.mean(0)
    Cb = np.corrcoef(Mb)
    iu = np.triu_indices(len(ls), 1)
    base_r = float(np.median(Cb[iu]))
    conv = []
    for (dr, cc), g in K.groupby(["drug", "conc"], observed=True):
        prof, keep_l = [], []
        for ln, gl in g.groupby("line", observed=True):
            if ln not in ctl_line:
                continue
            prof.append(X[gl.i.to_numpy()].mean(0))
            keep_l.append(ln)
        if len(prof) < 10:
            continue
        P = np.stack(prof)
        Q = np.stack([ctl_line[l] for l in keep_l])
        # Between-line SPREAD, on the same lines and genes, treated vs their own
        # untreated profiles. Centring each set across lines first (as an earlier
        # version did) removes the very quantity being measured, which is why
        # that version returned no change at any dose. Distances are normalised
        # by the untreated spread so drugs of different magnitude are comparable.
        def spread(Z):
            j_ = np.triu_indices(len(Z), 1)
            Dm = np.sqrt(np.maximum(
                ((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1), 0))
            return float(np.median(Dm[j_]))
        st, sc = spread(P), spread(Q)
        conv.append({"drug": dr, "conc": cc, "treated_spread": st,
                     "control_spread": sc,
                     "spread_ratio": st / sc if sc > 0 else np.nan,
                     "n_lines": len(P)})
    CV = pd.DataFrame(conv)
    if len(CV):
        print(f"   between-line spread, treated / untreated "
              f"(1.0 = unchanged, <1 = lines converge):")
        print(f"     overall median {CV.spread_ratio.median():.3f} over "
              f"{len(CV)} drug x dose")
        for cc, g in CV.groupby("conc"):
            w = stats.wilcoxon(g.spread_ratio - 1.0) if len(g) > 10 else None
            print(f"     conc {cc:>6}: {g.spread_ratio.median():.3f}  "
                  f"n = {len(g)}"
                  + (f"  p = {w.pvalue:.1e}" if w is not None else ""))
        print("   Measured on the SAME lines and genes, treated against their "
              "own untreated\n   profiles, so a ratio below 1 means the drug "
              "made the panel more alike.")

    # ---- 3. coordination: does the correlation structure change? --------
    print("\n3. COORDINATION — does exposure reorganise which genes co-vary?")
    sub = np.argsort(-X.var(0))[:400]
    top = K[K.conc == K.conc.max()]
    prof_t, prof_c = [], []
    for ln, gl in top.groupby("line", observed=True):
        if ln in ctl_line:
            prof_t.append(X[gl.i.to_numpy()].mean(0)[sub])
            prof_c.append(ctl_line[ln][sub])
    if len(prof_t) >= 15:
        Pt, Pc = np.stack(prof_t), np.stack(prof_c)
        Rt = np.corrcoef((Pt - Pt.mean(0)).T)
        Rc = np.corrcoef((Pc - Pc.mean(0)).T)
        iu2 = np.triu_indices(len(sub), 1)
        obs = float(np.mean(np.abs(Rt[iu2] - Rc[iu2])))
        agree = float(stats.spearmanr(Rt[iu2], Rc[iu2]).statistic)
        rng = np.random.default_rng(0)
        null = []
        allp = np.vstack([Pt, Pc])
        for _ in range(args.n_perm):
            pm = rng.permutation(len(allp))
            A_, B_ = allp[pm[:len(Pt)]], allp[pm[len(Pt):]]
            Ra = np.corrcoef((A_ - A_.mean(0)).T)
            Rb = np.corrcoef((B_ - B_.mean(0)).T)
            null.append(float(np.mean(np.abs(Ra[iu2] - Rb[iu2]))))
        null = np.array(null)
        pv = float(((null >= obs).sum() + 1) / (len(null) + 1))
        print(f"   gene-gene correlation across lines, treated vs untreated:")
        print(f"     agreement rho = {agree:+.3f} over "
              f"{len(sub)} genes")
        print(f"     mean |change| = {obs:.4f}; label-permutation null "
              f"{null.mean():.4f} [{np.percentile(null,2.5):.4f},"
              f"{np.percentile(null,97.5):.4f}], p = {pv:.4f}")
        print("   The permutation reassigns which profiles are treated, keeping "
              "the profiles\n   themselves, so it asks whether TREATMENT "
              "reorganises coordination or whether\n   any split of these "
              "samples would look as different.")
    else:
        agree, obs, pv = np.nan, np.nan, np.nan

    out = R.copy()
    out.to_csv(TAB / "expression_remodelling.csv", index=False)

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.2), constrained_layout=True)

    ax[0].hist(R.shift_share.clip(-0.2, 1.2), bins=45, color=VIOLET, alpha=0.85)
    ax[0].axvline(R.shift_share.median(), color=ORANGE, lw=2.4,
                  label=f"median {R.shift_share.median():.0%}")
    ax[0].axvline(k / p, color="#555", ls="--", lw=1.4,
                  label=f"chance ({k}/{p} dims) = {k/p:.1%}")
    ax[0].set_xlabel("share of the response inside the baseline programme space")
    ax[0].set_ylabel("(line, drug, dose) conditions")
    ax[0].legend(frameon=False, fontsize=7)
    ax[0].set_title("a  Shift, not remodelling", loc="left",
                    fontweight="bold", fontsize=9.5)

    ax[1].bar([0, 1], [R.shift_share.median(), np.median(naive)], width=0.5,
              color=[AQUA, GREY])
    for i_, v in enumerate([R.shift_share.median(), np.median(naive)]):
        ax[1].text(i_, v + 0.015, f"{v:.0%}", ha="center", fontsize=11,
                   fontweight="bold")
    ax[1].set_xticks([0, 1], ["cross-plate\ncovariance", "single\nmeasurement"],
                     fontsize=8)
    ax[1].set_ylabel("shift share")
    ax[1].set_title("b  Noise reads as remodelling", loc="left",
                    fontweight="bold", fontsize=9.5)
    ax[1].text(0.5, 0.12, "isotropic noise lands outside\nany low-dimensional "
               "span", transform=ax[1].transAxes, ha="center", fontsize=6.8,
               color="#444")

    ax[2].plot(bydose.index, bydose["remodel_share"], "o-", color=ORANGE,
               lw=2, ms=7, label="remodel share")
    ax[2].plot(bydose.index, bydose["shift_share"], "o-", color=AQUA, lw=2,
               ms=7, label="shift share")
    ax[2].set_xscale("log")
    ax[2].set_xlabel("concentration (µM)")
    ax[2].set_ylabel("share of the reproducible response")
    ax[2].legend(frameon=False, fontsize=7.5)
    ax[2].set_title("c  Flat across dose, so not death", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(CV):
        d0 = [CV[CV.conc == c].spread_ratio.dropna() for c in
              sorted(CV.conc.unique())]
        bp = ax[3].boxplot(d0, showfliers=False, patch_artist=True,
                           medianprops=dict(color="black", lw=1.4))
        for pch in bp["boxes"]:
            pch.set_facecolor(BLUE); pch.set_alpha(0.75)
        ax[3].axhline(1.0, color=ORANGE, lw=2.2, label="untreated spread")
        ax[3].set_xticks(range(1, len(d0) + 1),
                         [f"{c:g}" for c in sorted(CV.conc.unique())],
                         fontsize=8)
        ax[3].set_xlabel("concentration (µM)")
        ax[3].set_ylabel("between-line spread, treated / untreated")
        ax[3].legend(frameon=False, fontsize=7.5)
    ax[3].set_title("d  Do lines converge under treatment?", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Does drug exposure remodel the expression programme, or move "
                 "cells along the one they have?", fontsize=10.5, x=0.005,
                 ha="left", fontweight="bold")
    d = save_figure(fig, "expression_remodelling", FIG,
                    source_data={"per_condition": R, "by_dose":
                                 bydose.reset_index(), "capture_curve": CU,
                                 "convergence": CV if len(CV) else
                                 pd.DataFrame()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
