#!/usr/bin/env python3
"""What actually changes in a cell's transcriptome when a drug is added?

Section 30 established the geometry: exposure moves cells into directions that
untreated cells do not occupy, without converging cell identity or reorganising
co-expression. It never said what those directions ARE. This script names them.

The response of one condition is split three ways, mirroring the viability
decomposition but in gene space:

    delta(line, drug, dose) = common(dose)          every drug, every line
                            + signature(drug, dose) this drug, every line
                            + interaction(line,drug) this pairing only
                            + noise

and the question "what changed" has a different answer at each level:

  COMMON      the programme a cell runs because it was perturbed at all. If this
              dominates, then "the drug changed expression" is mostly a stress
              response and says little about the drug.
  SIGNATURE   the part that identifies the compound. This is what a mechanism of
              action should look like in transcription, and it is testable: two
              drugs with the same annotated mechanism should have more similar
              signatures than two drugs picked at random.
  INTERACTION the part specific to a (line, drug) pairing -- the quantity the
              rest of this project measures, here resolved into genes.

**Every term is measured as a covariance between plate 6 and plate 14**, Tahoe's
designed replicate pair. A variance would count noise as signal, and in gene
space that matters more than anywhere else: with 62,710 genes, noise has far more
room to look like structure than any real programme does. The same device makes
the per-gene numbers interpretable -- a gene's "response" is the part that
reproduces on an independent plate, not the part that moved once.

Three questions, each with the control that makes it answerable:

  1. How much of the response is common, signature and interaction? Reported per
     gene as well as globally, so "which genes are drug-discriminative" is
     answerable rather than asserted.
  2. What is the common programme biologically? Hallmark / Reactome / GO
     enrichment against a background of expressed genes, not all genes -- the
     usual error, which makes every result look like "expressed in cells".
  3. Do annotated mechanisms of action have distinguishable signatures? Within-
     MOA signature similarity against a label-permutation null, which is the
     only way to tell a real mechanism grouping from the fact that similar
     compounds were screened together.

Outputs: results/tables/what_changed_genes.csv
         results/tables/what_changed_moa.csv
         figure bundle results/figures/00_manuscript/what_changed/
"""
import argparse
import re
from collections import defaultdict
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
GS = ROOT / "data" / "external" / "genesets"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
CTRL = "DMSO_TF"
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
    # Keep genes expressed well enough that a fold change is meaningful. This
    # set is also the enrichment BACKGROUND: testing against all 62,710 genes
    # would report "expressed in cells" for every programme.
    expressed = X.mean(0) > 0.05
    print(f"{X.shape[0]} pseudobulks; {int(expressed.sum())} expressed genes "
          f"of {X.shape[1]}", flush=True)
    return C, X[:, expressed], G[expressed].reset_index(drop=True)


def read_gmt(path, universe, min_n=8, max_n=500):
    sets = {}
    if not path.exists():
        return sets
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 3:
            continue
        g = {x.strip().upper() for x in f[2:] if x.strip()} & universe
        if min_n <= len(g) <= max_n:
            sets[f[0]] = g
    return sets


def enrich(ranked, sets, universe, top_n=300):
    """Hypergeometric enrichment of the top-|loading| genes.

    The background is the set of EXPRESSED genes, not the genome; against the
    genome every programme in a cell line looks enriched for being expressed.
    """
    hit = set(ranked[:top_n])
    N, K = len(universe), len(hit)
    rows = []
    for name, g in sets.items():
        k = len(hit & g)
        if k < 3:
            continue
        p = stats.hypergeom.sf(k - 1, N, len(g), K)
        rows.append({"set": name, "n_set": len(g), "n_hit": k,
                     "expected": len(g) * K / N, "p": p})
    if not rows:
        return pd.DataFrame()
    E = pd.DataFrame(rows).sort_values("p")
    E["q"] = np.minimum(E.p * len(E) / np.arange(1, len(E) + 1).astype(float), 1)
    E["q"] = E.q[::-1].cummin()[::-1]
    return E


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--top-n", type=int, default=300)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    C, X, G = load()
    sym = G.gene_symbol.astype(str).str.upper().to_numpy()

    # plate-matched control, as the atlas's own pipeline does
    ctl = {}
    for (ln, pl), g in C[C.drug.astype(str) == CTRL].groupby(
            ["cell_line_id", "plate"], observed=True):
        ctl[(ln, pl)] = X[g.index.to_numpy()].mean(0)
    trt = C[C.drug.astype(str) != CTRL]
    rows = [(r.cell_line_id, r.drug, r.conc, r.plate, i)
            for i, r in zip(trt.index.to_numpy(), trt.itertuples())
            if (r.cell_line_id, r.plate) in ctl]
    K = pd.DataFrame(rows, columns=["line", "drug", "conc", "plate", "i"])
    D = np.stack([X[r.i] - ctl[(r.line, r.plate)] for r in K.itertuples()])
    print(f"{len(K)} treated conditions, {K.drug.nunique()} drugs, "
          f"{K.line.nunique()} lines", flush=True)

    # ---------- 1. three-level split, from cross-plate covariance ----------
    # Each level is an inner product between the SAME quantity measured on two
    # independent plates, so noise contributes zero in expectation to all three.
    rep = []
    for (ln, dr, cc), g in K.groupby(["line", "drug", "conc"], observed=True):
        if g.plate.nunique() < 2:
            continue
        pls = list(g.plate.unique())[:2]
        a = D[g.index[g.plate == pls[0]][0]]
        b = D[g.index[g.plate == pls[1]][0]]
        rep.append((ln, dr, cc, a, b))
    print(f"{len(rep)} conditions on two plates", flush=True)
    if len(rep) < 50:
        print("too few replicated conditions"); return

    A = np.stack([r[3] for r in rep])
    B = np.stack([r[4] for r in rep])
    meta = pd.DataFrame([(r[0], r[1], r[2]) for r in rep],
                        columns=["line", "drug", "conc"])

    # common programme: the mean response, measured on each plate separately
    cA, cB = A.mean(0), B.mean(0)
    # Drug signature per (drug, DOSE), not per drug. The same compound at 0.05
    # and 5 uM produces different responses in every line alike, and that
    # difference is a drug effect, not an interaction. Pooling doses leaves it
    # in the residual, where the line-effect term cannot reach it either, and it
    # is then reported as context-dependence -- which is what made a first
    # version of this script disagree with the scalar analysis of sec.31 by
    # sixty-fold on the same atlas.
    sigA, sigB = {}, {}
    for (dr, cc), g in meta.groupby(["drug", "conc"], observed=True):
        ii = g.index.to_numpy()
        if len(ii) < 2:
            continue
        sigA[(dr, cc)] = A[ii].mean(0) - cA
        sigB[(dr, cc)] = B[ii].mean(0) - cB
    keys = sorted(sigA)
    # per-DRUG signature, for the mechanism comparison: average its doses
    bydrug_A, bydrug_B = {}, {}
    for (dr, cc) in keys:
        bydrug_A.setdefault(dr, []).append(sigA[(dr, cc)])
        bydrug_B.setdefault(dr, []).append(sigB[(dr, cc)])
    drugs = sorted(bydrug_A)
    SA = np.stack([np.mean(bydrug_A[d], axis=0) for d in drugs])
    SB = np.stack([np.mean(bydrug_B[d], axis=0) for d in drugs])

    # per-gene variance at each level, each as a plate-by-plate covariance
    v_common = cA * cB
    v_sig = np.mean([sigA[k_] * sigB[k_] for k_ in keys], axis=0)

    # FOUR levels, not three. Subtracting the common programme and the drug
    # signature leaves the line's own general response sitting in the residual,
    # where it would be reported as interaction -- the exact error this project
    # corrected in RESULTS.md sec.27, and the first version of this script made
    # it too, returning 37% interaction against the 0.5% the scalar analysis
    # gives on the same atlas. The line term is estimated leave-one-DRUG-out so
    # a real interaction is not absorbed, and separately on each plate so the
    # two carry independent noise.
    dk = list(zip(meta.drug, meta.conc))
    r1A = A - cA - np.stack([sigA.get(k_, 0 * cA) for k_ in dk])
    r1B = B - cB - np.stack([sigB.get(k_, 0 * cB) for k_ in dk])

    def line_effect(R):
        out = np.zeros_like(R)
        for ln, g in meta.groupby("line", observed=True):
            ii = g.index.to_numpy()
            dv = meta.drug.to_numpy()[ii]
            by = {}
            for d in np.unique(dv):
                by[d] = R[ii[dv == d]].mean(0)
            if len(by) < 3:
                continue
            tot_l = np.sum(list(by.values()), axis=0)
            n_l = len(by)
            for i, d in zip(ii, dv):
                out[i] = (tot_l - by[d]) / (n_l - 1)
        return out

    aA, aB = line_effect(r1A), line_effect(r1B)
    v_line = (aA * aB).mean(0)
    gA, gB = r1A - aA, r1B - aB
    v_int = (gA * gB).mean(0)

    parts = [v_common, v_sig, v_line, v_int]
    tot = sum(np.maximum(v, 0) for v in parts)
    ok = tot > 0
    f_common, f_sig, f_line, f_int = [
        float(np.sum(np.maximum(v, 0)[ok]) / np.sum(tot[ok])) for v in parts]
    print(f"\n1. WHAT LEVEL DOES THE RESPONSE LIVE AT? "
          f"(reproducible variance, {len(drugs)} drugs)")
    print(f"   common programme (any drug, any line)   {f_common:6.1%}")
    print(f"   drug signature   (this drug, any line)  {f_sig:6.1%}")
    print(f"   cell property    (this line, any drug)  {f_line:6.1%}")
    print(f"   interaction      (this line x drug)     {f_int:6.1%}")
    print("   Each term is a plate-6 x plate-14 inner product, so noise "
          "contributes zero.")

    # Coherence check against sec.31, which put the same atlas's interaction at
    # 0.5% [0.0-1.5%]. Two numbers that far apart for a nominally similar
    # quantity have to be explained rather than reported side by side, so the
    # candidate explanation -- gene set -- is TESTED here rather than asserted.
    # Section 31 works on the 2,000 most response-variable genes; those are
    # where the shared drug effect is concentrated, which inflates its
    # denominator relative to a broad gene set.
    shared_share = f_int / max(f_common + f_sig + f_int, 1e-9)
    resp_idx = np.argsort(-A.var(0))[:2000]
    def split_on(idx):
        pr = [v[idx] for v in parts]
        t = sum(np.maximum(v, 0) for v in pr)
        m = t > 0
        f = [float(np.sum(np.maximum(v, 0)[m]) / np.sum(t[m])) for v in pr]
        return f[3] / max(f[0] + f[1] + f[3], 1e-9), f
    sh_resp, f_resp = split_on(resp_idx)
    rest = np.setdiff1d(np.arange(len(sym)), resp_idx)
    sh_rest, _ = split_on(rest)
    print(f"\n   coherence with sec.31 (which reports 0.5% [0.0-1.5%] on this "
          f"atlas):")
    print(f"     this split, all {len(sym)} expressed genes      "
          f"{shared_share:5.1%}   (shared+interaction denominator)")
    print(f"     restricted to the 2,000 response-variable genes  "
          f"{sh_resp:5.1%}")
    print(f"     the remaining {len(rest)} genes                  "
          f"{sh_rest:5.1%}")
    if sh_resp < 0.5 * shared_share:
        print("     The gene set explains much of the gap: the interaction is "
              "concentrated\n     OUTSIDE the response-variable genes that "
              "sec.31 selects, so a scalar over\n     those genes sees little "
              "of it. Both numbers are right about their own\n     question.")
    else:
        print("     The gene set does NOT explain the gap. The two analyses "
              "disagree on the\n     same atlas and the discrepancy is "
              "unresolved; neither number should be\n     quoted as the "
              "atlas's interaction share until it is.")

    # ---------- 2. what IS the common programme? ----------
    com = 0.5 * (cA + cB)
    order = np.argsort(-np.abs(com))
    up = [sym[i] for i in order if com[i] > 0][:25]
    dn = [sym[i] for i in order if com[i] < 0][:25]
    print(f"\n2. THE COMMON PROGRAMME — what every drug does to every line")
    print(f"   up:   {', '.join(up[:14])}")
    print(f"   down: {', '.join(dn[:14])}")
    rep_r = float(np.corrcoef(cA, cB)[0, 1])
    print(f"   reproducibility across plates: r = {rep_r:.3f}")

    universe = set(sym)
    ranked_up = [sym[i] for i in order if com[i] > 0]
    ranked_dn = [sym[i] for i in order if com[i] < 0]
    enr = {}
    for lib, fn in (("Hallmark", "MSigDB_Hallmark_2020.gmt"),
                    ("Reactome", "Reactome_2022.gmt"),
                    ("GO_BP", "GO_Biological_Process_2023.gmt")):
        sets = read_gmt(GS / fn, universe)
        if not sets:
            continue
        for tag, ranked in (("up", ranked_up), ("down", ranked_dn)):
            e = enrich(ranked, sets, universe, args.top_n)
            if len(e):
                e = e.assign(library=lib, direction=tag)
                enr[(lib, tag)] = e
                sig = e[e.q < 0.05]
                print(f"   {lib:9s} {tag:4s}: {len(sig)} sets at q<0.05"
                      + (f" — {'; '.join(sig.set.head(3))}" if len(sig) else ""))
    ENR = pd.concat(enr.values(), ignore_index=True) if enr else pd.DataFrame()

    # is the common programme the off-baseline axis of RESULTS sec.30?
    ctl_lines = sorted({k_[0] for k_ in ctl})
    M = np.stack([np.mean([ctl[k_] for k_ in ctl if k_[0] == l], axis=0)
                  for l in ctl_lines])
    Mc = M - M.mean(0, keepdims=True)
    _, _, Vt = np.linalg.svd(Mc, full_matrices=False)
    for k_ in (10, 20):
        Vk = Vt[:min(k_, Vt.shape[0])]
        inside = float(((Vk @ cA) @ (Vk @ cB)) / (cA @ cB))
        print(f"   share of the common programme inside the top-{k_} baseline "
              f"directions: {inside:.1%}")
    print("   A low share means the common programme is the new axis: the "
          "direction cells\n   move under treatment is not one on which "
          "untreated lines differ.")

    # ---------- 3. do mechanisms have distinguishable signatures? ----------
    md = pd.read_csv(TAB / "drug_metadata.csv")
    moa = dict(zip(md.drug.astype(str), md["moa-fine"].astype(str)))
    lab = np.array([moa.get(d, "nan") for d in drugs])
    good = ~np.isin(lab, ["nan", "unclear", "", "None"])
    # signature reproducibility per drug, the precondition for any MOA claim
    rep_drug = np.array([float(np.corrcoef(SA[i], SB[i])[0, 1])
                         for i in range(len(drugs))])
    print(f"\n3. DRUG SIGNATURES — do annotated mechanisms group?")
    print(f"   per-drug signature reproducibility (plate 6 vs 14): "
          f"median r = {np.median(rep_drug):+.3f}, "
          f"{(rep_drug > 0).mean():.0%} positive")

    # cross-plate similarity: drug i on plate A vs drug j on plate B. Using
    # different plates for the two members means a shared plate artefact cannot
    # make two drugs look alike.
    Za = SA / np.maximum(np.linalg.norm(SA, axis=1, keepdims=True), 1e-9)
    Zb = SB / np.maximum(np.linalg.norm(SB, axis=1, keepdims=True), 1e-9)
    Sim = Za @ Zb.T
    np.fill_diagonal(Sim, np.nan)
    idx = np.where(good)[0]
    counts = pd.Series(lab[idx]).value_counts()
    keep_moa = set(counts[counts >= 3].index)
    sel = np.array([i for i in idx if lab[i] in keep_moa])
    print(f"   {len(sel)} drugs across {len(keep_moa)} mechanisms with >=3 "
          f"members")
    if len(sel) > 20:
        L = lab[sel]
        S = Sim[np.ix_(sel, sel)]
        same = np.equal.outer(L, L) & ~np.eye(len(sel), dtype=bool)
        diff = ~np.equal.outer(L, L)
        obs = float(np.nanmean(S[same]) - np.nanmean(S[diff]))
        rng = np.random.default_rng(0)
        null = []
        for _ in range(args.n_perm):
            Lp = rng.permutation(L)
            sm = np.equal.outer(Lp, Lp) & ~np.eye(len(sel), dtype=bool)
            df = ~np.equal.outer(Lp, Lp)
            null.append(float(np.nanmean(S[sm]) - np.nanmean(S[df])))
        null = np.array(null)
        pv = float(((null >= obs).sum() + 1) / (len(null) + 1))
        print(f"   within-MOA minus between-MOA signature similarity: "
              f"{obs:+.4f}")
        print(f"   label-permutation null {null.mean():+.4f} "
              f"[{np.percentile(null,2.5):+.4f},{np.percentile(null,97.5):+.4f}]"
              f"   p = {pv:.4f}")
        # per-MOA, so the global number is not carried by one class
        per = []
        for m in sorted(keep_moa):
            w = np.where(L == m)[0]
            if len(w) < 3:
                continue
            sm = np.nanmean(S[np.ix_(w, w)][~np.eye(len(w), dtype=bool)])
            bt = np.nanmean(S[np.ix_(w, np.setdiff1d(np.arange(len(sel)), w))])
            per.append({"moa": m, "n_drugs": len(w), "within": float(sm),
                        "between": float(bt), "excess": float(sm - bt)})
        P = pd.DataFrame(per).sort_values("excess", ascending=False)
        print(f"\n   most transcriptionally distinct mechanisms:")
        print(P.head(8).round(4).to_string(index=False))
        P.to_csv(TAB / "what_changed_moa.csv", index=False)
    else:
        P, obs, pv, null = pd.DataFrame(), np.nan, np.nan, np.array([0.0])

    # ---------- per-gene table ----------
    GT = pd.DataFrame({
        "gene": sym, "common": com, "v_common": v_common, "v_signature": v_sig,
        "v_line": v_line, "v_interaction": v_int, "v_total": tot,
        "frac_signature": np.where(tot > 0, np.maximum(v_sig, 0) / np.maximum(
            tot, 1e-12), np.nan)})
    GT = GT.sort_values("v_common", ascending=False)
    GT.to_csv(TAB / "what_changed_genes.csv", index=False)
    # A fraction is only interpretable where the denominator is real signal:
    # ranked without this cut the list fills with pseudogenes whose total
    # reproducible variance is near zero and whose ratio is therefore arbitrary.
    cut = float(np.quantile(tot[ok], 0.90))
    disc = GT[(GT.v_signature > 0) & (GT.v_total >= cut)].sort_values(
        "frac_signature", ascending=False)
    print(f"\n   genes whose response is most DRUG-SPECIFIC rather than common:")
    print(f"   {', '.join(disc.gene.head(18))}")

    # ---------- figure ----------
    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 4, figsize=(18.5, 4.3), constrained_layout=True)

    vals = [f_common, f_sig, f_line, f_int]
    ax[0].bar(range(4), vals, width=0.6, color=[VIOLET, BLUE, GREY, ORANGE])
    for i_, v in enumerate(vals):
        ax[0].text(i_, v + 0.012, f"{v:.1%}", ha="center", fontsize=10,
                   fontweight="bold")
    ax[0].set_xticks(range(4), ["common\n(any drug)", "drug\nsignature",
                                "cell property\n(any drug)",
                                "line × drug\ninteraction"], fontsize=7)
    ax[0].set_ylabel("share of reproducible response variance")
    ax[0].set_ylim(0, max(vals) * 1.25)
    ax[0].set_title("a  Where the response lives", loc="left",
                    fontweight="bold", fontsize=9.5)

    n_show = 12
    o2 = np.argsort(-np.abs(com))
    gsel = [i for i in o2[:80]]
    gsel = sorted(gsel, key=lambda i: com[i])
    gsel = gsel[:n_show] + gsel[-n_show:]
    yy = np.arange(len(gsel))
    ax[1].barh(yy, [com[i] for i in gsel],
               color=[ORANGE if com[i] > 0 else AQUA for i in gsel], height=0.72)
    ax[1].set_yticks(yy, [sym[i] for i in gsel], fontsize=6.2)
    ax[1].axvline(0, color="#444", lw=0.9)
    ax[1].set_xlabel("common response (log₁₀ fold change)")
    ax[1].set_title(f"b  The programme every drug runs\n    (r = {rep_r:.2f} "
                    f"across plates)", loc="left", fontweight="bold",
                    fontsize=9.5)

    if len(ENR):
        e = ENR[(ENR.library == "Hallmark") & (ENR.q < 0.05)].head(10)
        if not len(e):
            e = ENR[ENR.q < 0.05].head(10)
        if len(e):
            yy2 = np.arange(len(e))[::-1]
            ax[2].barh(yy2, -np.log10(np.maximum(e.p, 1e-300)),
                       color=[ORANGE if d == "up" else AQUA
                              for d in e.direction], height=0.7)
            ax[2].set_yticks(yy2, [s[:38] for s in e.set], fontsize=6.0)
            ax[2].set_xlabel("−log₁₀ p (hypergeometric)")
    ax[2].set_title("c  What that programme is", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(P):
        top = P.head(12)
        yy3 = np.arange(len(top))[::-1]
        ax[3].barh(yy3, top.excess, color=VIOLET, height=0.7)
        ax[3].set_yticks(yy3, [f"{m[:30]} ({n})" for m, n in
                               zip(top.moa, top.n_drugs)], fontsize=6.0)
        ax[3].axvline(0, color="#444", lw=0.9)
        ax[3].set_xlabel("within-MOA − between-MOA signature similarity")
        ax[3].text(0.97, 0.05, f"overall {obs:+.3f}\npermutation p = {pv:.3f}",
                   transform=ax[3].transAxes, ha="right", fontsize=7,
                   color="#444")
    ax[3].set_title("d  Mechanisms with distinct signatures", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("What changes in a cell's transcriptome when a drug is added, "
                 "and how much of it is about the drug", fontsize=10.5,
                 x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "what_changed", FIG,
                    source_data={"per_gene": GT.head(4000),
                                 "enrichment": ENR, "per_moa": P},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
