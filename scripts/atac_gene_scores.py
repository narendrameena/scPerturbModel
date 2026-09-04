#!/usr/bin/env python3
"""The chromatin arm, retried on gene scores instead of motif deviations.

RESULTS.md sec.35 attempted the decomposition on Spear-ATAC ChromVar motif
deviations and could not use it: across all 123 (line, perturbation) pairs, zero
rejected the E-test (Peidli et al., *Nature Methods* 2024), so the perturbations
produced no detectable population shift and no interaction could be measured
either way. ChromVar deviations are a heavily compressed, already-regressed
summary -- 2,174 motifs standing in for the whole genome -- so the natural
question is whether the signal is there in a less lossy representation.

**Gene scores** are ArchR's per-cell, per-gene accessibility summary: 24,919
genes rather than 2,174 motifs, derived from the reads themselves rather than
from motif enrichment. If a single transcription-factor knockout moves chromatin
at all, this is where it should be visible.

The two questions from sec.35 are asked again in that space, unchanged:

  1. Do the perturbations move the cells? E-distance and the permutation E-test,
     per (line, perturbation), Bonferroni-corrected. This decides whether
     question 2 is answerable at all.
  2. If so, how does the response split into a common programme, a perturbation
     signature, a cell property and an interaction -- by the replicate-covariance
     estimator and by ANOVA variance components (Searle et al. 1992), the
     balanced-design limit of variancePartition (Hoffman & Schadt 2016).

**Why this is streamed.** The three matrices carry 589 million non-zero entries
between them, 36% dense for K562. Materialising them costs tens of gigabytes and
scipy's MatrixMarket reader is a text parser. Instead the file is read in chunks
and group sums are accumulated directly, so peak memory is the size of the
pseudobulk (a few hundred groups x 24,919 genes) rather than the size of the
data. A subsample of individual cells is retained in the same pass for the
E-test, which needs cells rather than pseudobulk.

Outputs: results/tables/atac_gene_scores.csv
         results/tables/atac_gene_scores_etest.csv
         figure bundle results/figures/00_manuscript/atac_gene_scores/
"""
import argparse
import gzip
import io
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
AT = ROOT / "data" / "external" / "scperturb_atac"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
LINES = ("K562", "GM12878", "MCF7")
CTRL_TOKENS = ("NTC", "CONTROL", "SAFE", "SCRAMBLE", "NEG")
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def stream_line(line, min_cells=15, n_cell_sample=2500, chunk=20_000_000,
                seed=0):
    """Pseudobulk one line by (perturbation, sample), plus a cell subsample.

    One pass over the MatrixMarket file. Group sums are accumulated with
    ``np.bincount`` on a flattened (group, gene) index, which is a C loop rather
    than the Python-level scatter ``np.add.at`` would use. Nothing larger than
    the pseudobulk is ever held.
    """
    z = zipfile.ZipFile(AT / f"PierceGreenleaf2021_{line}.zip")
    root = f"PierceGreenleaf2021_{line}/gene_scores"
    with z.open(f"{root}/obs.csv") as fh:
        obs = pd.read_csv(io.TextIOWrapper(fh), low_memory=False)
    with z.open(f"{root}/var.csv") as fh:
        var = pd.read_csv(io.TextIOWrapper(fh), low_memory=False)
    genes = var.iloc[:, 0].astype(str).to_numpy()
    obs["pert"] = obs.perturbation.astype(str)
    obs["rep"] = obs["sample"].astype(str)

    keys = obs.groupby(["pert", "rep"], observed=True).size()
    keys = keys[keys >= min_cells]
    gmap = {k: i for i, k in enumerate(keys.index)}
    gid = np.full(len(obs), -1, dtype=np.int64)
    for i, (p, r) in enumerate(zip(obs.pert, obs.rep)):
        gid[i] = gmap.get((p, r), -1)
    nG, nGene = len(gmap), len(genes)

    # cells kept for the E-test: controls plus a stratified sample of the rest,
    # so every perturbation contributes cells and controls are not swamped
    rng = np.random.default_rng(seed)
    isc = obs.pert.str.upper().str.contains("|".join(CTRL_TOKENS), na=False)
    pick = list(rng.choice(np.where(isc)[0],
                           min(int(isc.sum()), n_cell_sample), replace=False))
    per_p = max(n_cell_sample // max(obs.pert.nunique(), 1), 30)
    for p, g in obs[~isc].groupby("pert", observed=True):
        ii = g.index.to_numpy()
        pick += list(rng.choice(ii, min(len(ii), per_p), replace=False))
    pick = np.array(sorted(set(pick)))
    prow = np.full(len(obs), -1, dtype=np.int64)
    prow[pick] = np.arange(len(pick))

    acc = np.zeros(nG * nGene, dtype=np.float64)
    cnt = np.zeros(nG, dtype=np.int64)
    for i in range(len(obs)):
        if gid[i] >= 0:
            cnt[gid[i]] += 1
    sub = np.zeros((len(pick), nGene), dtype=np.float32)

    nnz = 0
    with z.open(f"{root}/counts.mtx.gz") as fh:
        with gzip.open(fh, "rt") as gz:
            header = 0
            for ln_ in gz:
                if ln_.startswith("%"):
                    continue
                header += 1
                break
            for ch in pd.read_csv(gz, sep=r"\s+", header=None,
                                  names=["r", "c", "v"], chunksize=chunk,
                                  engine="c", dtype={"r": np.int64,
                                                     "c": np.int64,
                                                     "v": np.float64}):
                r = ch.r.to_numpy() - 1
                c = ch.c.to_numpy() - 1
                v = ch.v.to_numpy()
                nnz += len(r)
                g = gid[r]
                m = g >= 0
                if m.any():
                    acc += np.bincount(g[m] * nGene + c[m], weights=v[m],
                                       minlength=nG * nGene)
                pr = prow[r]
                m2 = pr >= 0
                if m2.any():
                    sub[pr[m2], c[m2]] = v[m2]
    PB = (acc.reshape(nG, nGene) /
          np.maximum(cnt, 1)[:, None]).astype(np.float32)
    K = pd.DataFrame(list(gmap.keys()), columns=["pert", "rep"])
    print(f"  {line:8s} {len(obs):7,} cells x {nGene:,} genes, {nnz:,} nnz -> "
          f"{nG} pseudobulks, {len(pick):,} cells kept for the E-test",
          flush=True)
    return K, PB, genes, sub, obs.iloc[pick].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--n-genes-etest", type=int, default=2000)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from published_methods_check import anova_components, etest
    from atac_decomposition import LINES as _L

    print("streaming Spear-ATAC gene scores ...", flush=True)
    PBs, cellsub, genes = {}, {}, None
    for ln in LINES:
        K, PB, g, sub, sobs = stream_line(ln)
        if genes is None:
            genes = g
        elif not np.array_equal(genes, g):
            raise SystemExit("gene sets differ between lines")
        PBs[ln] = (K, PB)
        cellsub[ln] = (sobs, sub)

    # control-subtracted deltas, control = non-targeting guides in the same
    # (line, sample), so a sample-level offset cancels
    rows = []
    for ln in LINES:
        K, PB = PBs[ln]
        isc = K.pert.str.upper().str.contains("|".join(CTRL_TOKENS), na=False)
        if not isc.any():
            print(f"  {ln}: no control guides; skipping"); continue
        ctl = {r: PB[g.index[isc[g.index]].to_numpy()].mean(0)
               for r, g in K.groupby("rep", observed=True)
               if isc[g.index].any()}
        for i in np.where(~isc)[0]:
            r = K.rep.iloc[i]
            if r in ctl:
                rows.append({"line": ln, "pert": K.pert.iloc[i], "rep": r,
                             "v": PB[i] - ctl[r]})
    R = pd.DataFrame(rows)
    shared = sorted(set.intersection(*[set(R[R.line == l].pert)
                                       for l in LINES]))
    R = R[R.pert.isin(shared)].reset_index(drop=True)
    V = np.stack(R.v.to_numpy())
    print(f"\n{len(R)} (line, perturbation, sample) rows, {len(shared)} shared "
          f"perturbations, {V.shape[1]:,} genes", flush=True)

    # ---------- 1. E-test first: is there anything to decompose? ----------
    print(f"\n1. E-distance and E-test (Peidli et al. 2024) on gene scores",
          flush=True)
    # rank genes on the CONTROL cells only, so the selection cannot favour any
    # perturbation and the test is not tuned on its own outcome
    ctl_all = []
    for ln in LINES:
        so, sx = cellsub[ln]
        m = so.pert.str.upper().str.contains("|".join(CTRL_TOKENS), na=False)
        if m.any():
            ctl_all.append(sx[m.to_numpy()])
    gsel = np.argsort(-np.vstack(ctl_all).var(0))[:args.n_genes_etest]
    et = []
    for ln in LINES:
        so, sx = cellsub[ln]
        isc = so.pert.str.upper().str.contains("|".join(CTRL_TOKENS), na=False)
        C = sx[np.where(isc)[0]][:, gsel]
        if len(C) < 20:
            continue
        for p, g in so[~isc].groupby("pert", observed=True):
            if p not in shared:
                continue
            Y = sx[g.index.to_numpy()][:, gsel]
            if len(Y) < 20:
                continue
            e, pv = etest(C, Y, n_perm=args.n_perm)
            et.append({"line": ln, "pert": p, "edist": e, "p": pv,
                       "n_cells": len(Y)})
    ET = pd.DataFrame(et)
    n_sig = 0
    if len(ET):
        ET["q"] = np.minimum(ET.p * len(ET), 1.0)
        n_sig = int((ET.q < 0.05).sum())
        print(f"   {len(ET)} (line, perturbation) tests on "
              f"{len(gsel)} genes; {n_sig} reject at Bonferroni q<0.05 "
              f"({n_sig/len(ET):.0%})")
        print(f"   median E-distance {ET.edist.median():.4f}, "
              f"max {ET.edist.max():.4f}")
        ET.to_csv(TAB / "atac_gene_scores_etest.csv", index=False)
        if n_sig == 0:
            print("   As with ChromVar motifs, no perturbation moves the "
                  "population. The\n   decomposition below is reported for "
                  "completeness but cannot be read as\n   evidence about "
                  "interactions.")
        else:
            print("   Unlike ChromVar motifs, perturbations ARE detectable "
                  "here, so the\n   decomposition below is interpretable.")

    # ---------- 2. the decomposition ----------
    lines_l = sorted(R.line.unique())
    perts_l = [p for p in shared
               if all(((R.line == l) & (R.pert == p)).sum() >= 2
                      for l in lines_l)]
    rng0 = np.random.default_rng(0)
    cube = np.empty((len(lines_l), len(perts_l), 2, V.shape[1]),
                    dtype=np.float64)
    for ia, l in enumerate(lines_l):
        for jb, p in enumerate(perts_l):
            ii = R.index[(R.line == l) & (R.pert == p)].to_numpy()
            cube[ia, jb] = V[rng0.choice(ii, 2, replace=False)]
    print(f"\n2. Variance components, ANOVA method of moments "
          f"(Searle et al. 1992)")
    print(f"   design {cube.shape[0]} lines x {cube.shape[1]} perturbations x "
          f"{cube.shape[2]} replicates x {cube.shape[3]:,} genes", flush=True)
    VP, neg = anova_components(cube, *cube.shape[:3])
    VP = VP.dropna()
    print(f"   context (cell line)         {VP.ctx.median():6.1%}")
    print(f"   perturbation                {VP.pert.median():6.1%}")
    print(f"   context x perturbation      {VP.inter.median():6.1%}")
    print(f"   residual (within-replicate) {VP.resid.median():6.1%}")
    print(f"   negative moment estimates: context {neg['ctx']:.0%}, "
          f"perturbation {neg['pert']:.0%}, interaction {neg['inter']:.0%}")
    rep = VP.ctx + VP.pert + VP.inter
    share = float(np.median(VP.inter[rep > 0] / rep[rep > 0]))
    print(f"   interaction as a share of reproducible variance: {share:.1%}")

    T = pd.DataFrame([{"assay": "gene_scores", "n_lines": len(lines_l),
                       "n_perturbations": len(perts_l),
                       "n_genes": cube.shape[3],
                       "ctx": VP.ctx.median(), "pert": VP.pert.median(),
                       "inter": VP.inter.median(), "resid": VP.resid.median(),
                       "interaction_of_reproducible": share,
                       "n_etest": len(ET), "n_etest_sig": n_sig,
                       "median_edist": float(ET.edist.median()) if len(ET)
                       else np.nan}])
    T.to_csv(TAB / "atac_gene_scores.csv", index=False)

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    if len(ET):
        ax[0].hist(ET.edist, bins=28, color=AQUA, alpha=0.85)
        ax[0].axvline(float(ET.edist.median()), color=ORANGE, lw=2.2,
                      label=f"median {ET.edist.median():.3f}")
        ax[0].set_xlabel("energy distance, perturbed vs control cells")
        ax[0].set_ylabel("(line, perturbation) pairs")
        ax[0].legend(frameon=False, fontsize=7.5)
        ax[0].text(0.5, 0.75, f"{n_sig} of {len(ET)} reject\n"
                   f"at Bonferroni q < 0.05", transform=ax[0].transAxes,
                   ha="center", fontsize=9.5, fontweight="bold",
                   color=ORANGE if n_sig == 0 else AQUA)
    ax[0].set_title("a  Do the perturbations move the cells?", loc="left",
                    fontweight="bold", fontsize=9.5)

    med = [VP.ctx.median(), VP.pert.median(), VP.inter.median(),
           VP.resid.median()]
    ax[1].bar(range(4), med, width=0.6, color=[GREY, BLUE, ORANGE, "#d8d8d8"])
    for i_, v in enumerate(med):
        ax[1].text(i_, v + 0.012, f"{v:.1%}", ha="center", fontsize=9,
                   fontweight="bold")
    ax[1].set_xticks(range(4), ["cell line", "perturbation",
                                "line ×\nperturbation", "residual"],
                     fontsize=7)
    ax[1].set_ylabel("variance fraction")
    ax[1].set_ylim(0, 1.05)
    ax[1].set_title("b  Variance components, gene scores", loc="left",
                    fontweight="bold", fontsize=9.5)

    prev = read_prev()
    xs = ["ChromVar\nmotifs (2,174)", f"gene scores\n({cube.shape[3]:,})"]
    ax[2].bar([0, 1], [prev, share], width=0.5, color=[GREY, VIOLET])
    for i_, v in enumerate([prev, share]):
        ax[2].text(i_, v + 0.012, f"{v:.1%}", ha="center", fontsize=10,
                   fontweight="bold")
    ax[2].set_xticks([0, 1], xs, fontsize=7.5)
    ax[2].set_ylabel("interaction / reproducible variance")
    ax[2].text(0.5, 0.85, "neither is interpretable unless\npanel a rejects",
               transform=ax[2].transAxes, ha="center", fontsize=7,
               color="#444")
    ax[2].set_title("c  Representation changes the number", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("The chromatin arm retried on gene scores rather than motif "
                 "deviations", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    d = save_figure(fig, "atac_gene_scores", FIG,
                    source_data={"summary": T, "etest": ET,
                                 "variance_components": VP}, script=__file__)
    print(f"figure bundle -> {d}")


def read_prev():
    f = TAB / "published_methods_check.csv"
    if f.exists():
        t = pd.read_csv(f)
        m = t[t.method.str.contains("variance", case=False)]
        if len(m):
            return float(m.interaction_of_reproducible.iloc[0])
    return 0.307


if __name__ == "__main__":
    main()
