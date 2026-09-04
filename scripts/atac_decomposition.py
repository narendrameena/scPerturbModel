#!/usr/bin/env python3
"""Does the decomposition hold in chromatin, and in a fully crossed design?

Everything measured so far has been transcription (Tahoe, LINCS, sci-Plex, OP3)
or viability (PRISM, GDSC). Both are downstream readouts of the same cell state,
so agreement between them is weaker evidence than it looks. Chromatin
accessibility is a different layer, and Spear-ATAC (Pierce, Greenleaf et al.
2021, distributed through scPerturb) supplies it in the one design that removes
most of the usual excuses:

  * **fully crossed** -- all 41 CRISPR perturbations are present in all three
    cell lines (K562, GM12878, MCF7), so no context x perturbation cell is
    missing and nothing has to be imputed or dropped;
  * **replicated** -- 4 to 6 independent samples per line, which is what makes
    the interaction estimable as a covariance rather than a variance;
  * **genetic, not chemical** -- so a result here is not a property of how
    compounds are dosed;
  * **motif space** -- ChromVar TF deviation scores are the same features in
    every cell line, so no gene-set intersection is needed and the three lines
    are directly comparable.

The same four-way split is applied as in gene space:

    y(line, perturbation) = common + signature(pert) + alpha(line)
                          + gamma(line, pert) + noise

with every term measured as a covariance between INDEPENDENT replicate samples,
so noise contributes zero to each. The two questions that matter:

  1. Is there a cell property in chromatin -- a line that responds to every
     perturbation? If the answer is no, the correction this project introduced
     is specific to viability and transcription rather than general.
  2. How large is the interaction once that property is removed? Tahoe says
     0.5% at matched dose in transcription; a fully crossed, replicated design
     in a different layer is the cleanest available check on whether that is a
     property of the biology or of that atlas.

Because only three contexts are available, alpha is estimated leave-one-
PERTURBATION-out and centred across the three lines, and the interaction is
reported with a permutation null over perturbation labels -- with 3 contexts a
point estimate alone would not be interpretable.

Outputs: results/tables/atac_decomposition.csv
         figure bundle results/figures/00_manuscript/atac_decomposition/
"""
import argparse
import io
import zipfile
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
AT = ROOT / "data" / "external" / "scperturb_atac"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
LINES = ("K562", "GM12878", "MCF7")
CTRL_TOKENS = ("NTC", "CONTROL", "SAFE", "SCRAMBLE", "NEG")
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def load_line(line, assay="ChromVar"):
    """Pseudobulk one cell line to (perturbation, sample) x feature.

    Cells are aggregated per (perturbation, sample) because the sample is the
    replicate unit: two cells carrying the same guide in the same well are not
    independent measurements of that perturbation's effect, and treating them as
    replicates would make noise look reproducible.
    """
    z = zipfile.ZipFile(AT / f"PierceGreenleaf2021_{line}.zip")
    root = f"PierceGreenleaf2021_{line}/{assay}"
    with z.open(f"{root}/obs.csv") as fh:
        obs = pd.read_csv(io.TextIOWrapper(fh), low_memory=False)
    with z.open(f"{root}/var.csv") as fh:
        var = pd.read_csv(io.TextIOWrapper(fh), low_memory=False)
    with z.open(f"{root}/counts.mtx.gz") as fh:
        import gzip
        M = mmread(io.BytesIO(gzip.decompress(fh.read())))
    # mmread returns a dense ndarray for "array"-format MatrixMarket files and
    # a sparse matrix for "coordinate" format; ChromVar deviations are dense.
    X = np.asarray(M.todense() if hasattr(M, "todense") else M,
                   dtype=np.float32)
    if X.shape[0] == len(var) and X.shape[1] == len(obs):
        X = X.T
    feat = var.iloc[:, 0].astype(str).to_numpy()
    obs["pert"] = obs.perturbation.astype(str)
    obs["rep"] = obs["sample"].astype(str)
    print(f"  {line:8s} {X.shape[0]:7,} cells x {X.shape[1]:5,} features, "
          f"{obs.pert.nunique()} perturbations, {obs.rep.nunique()} samples",
          flush=True)
    return obs, X, feat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assay", default="ChromVar")
    ap.add_argument("--min-cells", type=int, default=15)
    ap.add_argument("--n-perm", type=int, default=500)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    print(f"loading Spear-ATAC ({args.assay}) ...", flush=True)
    PB, feats = {}, None
    for ln in LINES:
        obs, X, feat = load_line(ln, args.assay)
        if feats is None:
            feats = feat
        elif not np.array_equal(feats, feat):
            common = [f for f in feats if f in set(feat)]
            idx_a = {f: i for i, f in enumerate(feats)}
            idx_b = {f: i for i, f in enumerate(feat)}
            X = X[:, [idx_b[f] for f in common]]
            for k in PB:
                PB[k] = (PB[k][0], PB[k][1][:, [idx_a[f] for f in common]])
            feats = np.array(common)
        keys, mats = [], []
        for (p, r), g in obs.groupby(["pert", "rep"], observed=True):
            ii = g.index.to_numpy()
            if len(ii) < args.min_cells:
                continue
            keys.append((p, r))
            mats.append(X[ii].mean(0))
        PB[ln] = (pd.DataFrame(keys, columns=["pert", "rep"]),
                  np.stack(mats).astype(np.float32))
        print(f"    -> {len(keys)} (perturbation, sample) pseudobulks",
              flush=True)

    # control: the non-targeting guides in the same line and sample, so a
    # sample-level batch offset cancels exactly as a plate-matched control does
    D = {}
    for ln in LINES:
        Kl, Xl = PB[ln]
        is_ctrl = Kl.pert.str.upper().str.contains("|".join(CTRL_TOKENS),
                                                   na=False)
        if not is_ctrl.any():
            print(f"  {ln}: no control guides found; using the per-sample mean "
                  f"over perturbations instead")
        ctl = {}
        for r, g in Kl.groupby("rep", observed=True):
            sel = g.index[is_ctrl[g.index]] if is_ctrl.any() else g.index
            ctl[r] = Xl[sel.to_numpy()].mean(0)
        keep = ~is_ctrl if is_ctrl.any() else np.ones(len(Kl), bool)
        rows, mats = [], []
        for i in np.where(keep)[0]:
            r = Kl.rep.iloc[i]
            rows.append((Kl.pert.iloc[i], r))
            mats.append(Xl[i] - ctl[r])
        D[ln] = (pd.DataFrame(rows, columns=["pert", "rep"]),
                 np.stack(mats).astype(np.float32))

    shared = set(D[LINES[0]][0].pert)
    for ln in LINES[1:]:
        shared &= set(D[ln][0].pert)
    shared = sorted(shared)
    print(f"\n{len(shared)} perturbations present in all {len(LINES)} lines, "
          f"{len(feats)} features", flush=True)
    if len(shared) < 10:
        print("too few shared perturbations"); return

    # two independent replicate halves per (line, perturbation)
    A, B, meta = [], [], []
    for ln in LINES:
        Kl, Xl = D[ln]
        for p, g in Kl.groupby("pert", observed=True):
            if p not in shared or g.rep.nunique() < 2:
                continue
            reps = sorted(g.rep.unique())
            h = max(len(reps) // 2, 1)
            ia = g.index[g.rep.isin(reps[:h])].to_numpy()
            ib = g.index[g.rep.isin(reps[h:])].to_numpy()
            if not len(ia) or not len(ib):
                continue
            A.append(Xl[ia].mean(0)); B.append(Xl[ib].mean(0))
            meta.append((ln, p))
    A, B = np.stack(A), np.stack(B)
    meta = pd.DataFrame(meta, columns=["line", "pert"])
    print(f"{len(meta)} (line, perturbation) conditions with two independent "
          f"replicate halves", flush=True)

    def split(Am, Bm, mt):
        cA, cB = Am.mean(0), Bm.mean(0)
        sA, sB = {}, {}
        for p, g in mt.groupby("pert", observed=True):
            ii = g.index.to_numpy()
            sA[p] = Am[ii].mean(0) - cA
            sB[p] = Bm[ii].mean(0) - cB
        r1A = Am - cA - np.stack([sA[p] for p in mt.pert])
        r1B = Bm - cB - np.stack([sB[p] for p in mt.pert])

        def line_eff(R):
            out = np.zeros_like(R)
            for ln, g in mt.groupby("line", observed=True):
                ii = g.index.to_numpy()
                pv = mt.pert.to_numpy()[ii]
                by = {p: R[ii[pv == p]].mean(0) for p in np.unique(pv)}
                if len(by) < 3:
                    continue
                tot, n = np.sum(list(by.values()), axis=0), len(by)
                for i, p in zip(ii, pv):
                    out[i] = (tot - by[p]) / (n - 1)
            return out
        aA, aB = line_eff(r1A), line_eff(r1B)
        gA, gB = r1A - aA, r1B - aB
        parts = [cA * cB,
                 np.mean([sA[p] * sB[p] for p in sA], axis=0),
                 (aA * aB).mean(0), (gA * gB).mean(0)]
        tot = sum(np.maximum(v, 0) for v in parts)
        ok = tot > 0
        return [float(np.sum(np.maximum(v, 0)[ok]) / np.sum(tot[ok]))
                for v in parts]

    f_com, f_sig, f_line, f_int = split(A, B, meta)
    print(f"\nFOUR-WAY SPLIT ({args.assay}, {len(LINES)} lines, "
          f"{len(shared)} perturbations)")
    print(f"  common programme   (any perturbation, any line)  {f_com:6.1%}")
    print(f"  perturbation sig.  (this perturbation, any line) {f_sig:6.1%}")
    print(f"  cell property      (this line, any perturbation) {f_line:6.1%}")
    print(f"  interaction        (this line x perturbation)    {f_int:6.1%}")

    # is the cell property real? Its estimate on disjoint halves of the
    # perturbations should agree, exactly the test used in PRISM (r = 0.989).
    perts = sorted(meta.pert.unique())
    rng = np.random.default_rng(0)
    pm = rng.permutation(len(perts))
    hA = {perts[i] for i in pm[:len(perts) // 2]}
    aa, bb = [], []
    for ln, g in meta.groupby("line", observed=True):
        ia = g.index[g.pert.isin(hA)].to_numpy()
        ib = g.index[~g.pert.isin(hA)].to_numpy()
        if len(ia) and len(ib):
            aa.append(0.5 * (A[ia] + B[ia]).mean(0))
            bb.append(0.5 * (A[ib] + B[ib]).mean(0))
    aa, bb = np.stack(aa), np.stack(bb)
    aa = aa - aa.mean(0); bb = bb - bb.mean(0)
    r_alpha = float(np.corrcoef(aa.ravel(), bb.ravel())[0, 1])
    print(f"\n  cell property reproduces across disjoint perturbation halves: "
          f"r = {r_alpha:+.3f} ({len(aa)} lines)")
    print("   With only three contexts this is a weak test -- it is reported "
          "because a\n   NEGATIVE result would matter, not because a positive "
          "one is conclusive.")

    # permutation null for the interaction: shuffle perturbation labels within
    # each line, which destroys any line x perturbation pairing while leaving
    # every marginal, the common programme and the cell property intact
    null = []
    for _ in range(args.n_perm):
        mt = meta.copy()
        idx = np.arange(len(mt))
        for ln, g in mt.groupby("line", observed=True):
            ii = g.index.to_numpy()
            idx[ii] = rng.permutation(ii)
        mt2 = meta.copy()
        mt2["pert"] = meta.pert.to_numpy()[idx]
        null.append(split(A, B, mt2)[3])
    null = np.array(null)
    pv = float(((null >= f_int).sum() + 1) / (len(null) + 1))
    print(f"\n  interaction {f_int:.1%} against a label-permutation null of "
          f"{null.mean():.1%} [{np.percentile(null,2.5):.1%},"
          f"{np.percentile(null,97.5):.1%}]   p = {pv:.4f}")

    T = pd.DataFrame([{"assay": args.assay, "n_lines": len(LINES),
                       "n_perturbations": len(shared), "n_features": len(feats),
                       "common": f_com, "signature": f_sig,
                       "cell_property": f_line, "interaction": f_int,
                       "alpha_split_r": r_alpha, "null_mean": float(null.mean()),
                       "p_vs_null": pv}])
    T.to_csv(TAB / "atac_decomposition.csv", index=False)

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    vals = [f_com, f_sig, f_line, f_int]
    ax[0].bar(range(4), vals, width=0.6, color=[VIOLET, BLUE, GREY, ORANGE])
    for i_, v in enumerate(vals):
        ax[0].text(i_, v + 0.012, f"{v:.1%}", ha="center", fontsize=10,
                   fontweight="bold")
    ax[0].set_xticks(range(4), ["common", "perturbation\nsignature",
                               "cell\nproperty", "line ×\nperturbation"],
                     fontsize=7)
    ax[0].set_ylabel("share of reproducible variance")
    ax[0].set_title(f"a  Chromatin, fully crossed\n    ({len(LINES)} lines × "
                    f"{len(shared)} perturbations)", loc="left",
                    fontweight="bold", fontsize=9.5)

    ax[1].hist(null, bins=30, color=GREY, alpha=0.85, label="label-permutation")
    ax[1].axvline(f_int, color=ORANGE, lw=2.4, label=f"observed {f_int:.1%}")
    ax[1].set_xlabel("interaction share")
    ax[1].set_ylabel("permutations")
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].set_title(f"b  Against its own null (p = {pv:.3f})", loc="left",
                    fontweight="bold", fontsize=9.5)

    ax[2].scatter(aa.ravel(), bb.ravel(), s=5, alpha=0.25, color=VIOLET,
                  edgecolors="none")
    ax[2].set_xlabel("cell property, perturbation half A")
    ax[2].set_ylabel("cell property, perturbation half B")
    ax[2].text(0.04, 0.94, f"r = {r_alpha:+.3f}\n{len(aa)} lines only",
               transform=ax[2].transAxes, va="top", fontsize=8,
               fontweight="bold")
    ax[2].set_title("c  Does a cell property exist in chromatin?", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("The decomposition in a third layer: chromatin accessibility, "
                 "genetic perturbation, fully crossed", fontsize=10.5, x=0.005,
                 ha="left", fontweight="bold")
    d = save_figure(fig, "atac_decomposition", FIG, source_data={"summary": T},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
