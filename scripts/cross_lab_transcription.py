#!/usr/bin/env python3
"""Does the cross-laboratory result hold for transcription, or only viability?

The viability arm (cross_lab_reproducibility.py) found that ~56% of the
reproducible line-specific drug response survives the move from PRISM (Broad,
pooled-barcode viability) to GDSC (Sanger, fitted IC50). That comparison changes
laboratory AND assay at once, and it measures a scalar phenotype. If the number
is a property of laboratories rather than of the killing assay, it should recur
for a transcriptional readout.

Two comparisons, matched in construction:

  within-lab   LINCS phase 1 (GSE92742) vs phase 2 (GSE70138) — both Broad, both
               L1000, different experiments years apart. 1,034 shared compounds,
               18 shared cell lines. This is the transcriptional ceiling.
  cross-lab    Tahoe-100M (Vevo/Arc, single-cell) vs LINCS (Broad, bulk L1000).
               ~133 shared compounds.

**The metric differs from the viability arm for a reason.** There the residual
was a scalar and could be correlated ACROSS cell lines, needing >=20 shared lines
per compound. Here only a handful of lines are shared, so correlating across
lines would be hopeless — but each (line, compound) carries a whole residual
GENE VECTOR, so the profile itself can be correlated over hundreds of genes. That
asks the same question ("does this line's compound-specific deviation reproduce")
with the power placed on the gene axis instead of the line axis.

Residuals are the project's standard construction: the line's profile minus the
leave-one-line-out mean profile for that compound, so the compound's shared
effect is removed and only the context-specific part is compared.

Outputs: results/tables/cross_lab_transcription.csv
         figure bundle results/figures/16_crosslab/cross_lab_transcription/
"""
import argparse
import gzip
import re
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.celldrug import remove_line_effect_profiles
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
LIN = ROOT / "data" / "external" / "lincs"
FIG = ROOT / "results" / "figures" / "16_crosslab"
TAB = ROOT / "results" / "tables"
MIN_GENES = 100
MIN_LINES_PER_CPD = 3
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def load_lincs(gctx, inst_info, gene_info, tag):
    """Per (cell line, compound) residual profile on the landmark genes."""
    gi = pd.read_csv(gene_info, sep="\t")
    lm = gi[gi.pr_is_lm == 1] if "pr_is_lm" in gi.columns else gi
    landmark = lm.pr_gene_id.astype(str).tolist()
    sym = dict(zip(gi.pr_gene_id.astype(str), gi.pr_gene_symbol.astype(str)))
    inst = pd.read_csv(inst_info, sep="\t", low_memory=False)
    cp = inst[inst.pert_type == "trt_cp"].copy()
    cp["k"] = cp.pert_iname.map(norm)
    keep = cp.groupby("k").cell_id.nunique().loc[lambda s: s >= 2].index
    cp = cp[cp.k.isin(keep)]
    print(f"  {tag}: {len(cp)} instances, {cp.k.nunique()} compounds, "
          f"{cp.cell_id.nunique()} lines", flush=True)
    from cmapPy.pandasGEXpress.parse_gctx import parse
    g = parse(str(gctx), rid=landmark, cid=cp.inst_id.tolist())
    M = g.data_df.T
    meta = cp.set_index("inst_id").loc[M.index]
    genes = np.array([sym.get(c, c) for c in M.columns])
    X = M.to_numpy(dtype=np.float32)
    out = {}
    for k, gg in meta.groupby("k", observed=True):
        pos = {c: i for i, c in enumerate(M.index)}
        by_line = {}
        acc = {}
        for cid, gl in gg.groupby("cell_id", observed=True):
            ii = [pos[x] for x in gl.index]
            # LINCS appends assay suffixes (A549.311, HELA.311). Stripping them
            # makes several ids collide -- ASC/ASC.C, NPC/NPC.CAS9/NPC.TAK,
            # SKL/SKL.C -- and plain dict assignment silently DISCARDED the
            # earlier groups, dropping 12,313 of 21,627 instances in colliding
            # groups. Accumulate and average instead. NPC.CAS9 is a Cas9
            # derivative rather than the same culture, so this is a compromise;
            # it is at least no longer a silent last-wins.
            k = norm(str(cid).split(".")[0])
            acc.setdefault(k, []).append(X[ii])
        for k, mats in acc.items():
            by_line[k] = np.concatenate(mats, axis=0).mean(0)
        if len(by_line) < 2:
            continue
        names = list(by_line)
        A = np.stack([by_line[c] for c in names])
        tot = A.sum(0)
        for j, c in enumerate(names):
            out[(c, k)] = A[j] - (tot - A[j]) / (len(names) - 1)
    # remove each line's general response before calling the rest a
    # cell-drug relation -- see celldrug.remove_line_effect_profiles
    return remove_line_effect_profiles(out), genes


def load_tahoe(pb_dir):
    from perturbmodel.evaluation.delta_eval import (build_deltas,
                                                    load_pseudobulk)
    X, cond = load_pseudobulk(ROOT / pb_dir)
    G, DELTA = build_deltas(X, cond, keep_plate=True)
    genes = np.array(cond["genes"]) if isinstance(cond, dict) and "genes" in cond \
        else None
    cl = pd.read_parquet(ROOT / "data/metadata/metadata/cell_line_metadata.parquet")
    cvcl2name = dict(zip(cl.Cell_ID_Cellosaur.astype(str),
                         cl.cell_name.astype(str)))
    # gene names live beside the counts, in the same column order as DELTA
    gnames = pd.read_csv(ROOT / pb_dir / "genes.csv")
    col = [c for c in ("gene_symbol", "symbol", "gene", "gene_name")
           if c in gnames.columns]
    genes = np.array(gnames[col[0]].astype(str)) if col else \
        np.array(gnames.iloc[:, 0].astype(str))
    K = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                      "drug": G.drug})
    K["cname"] = K.line.astype(str).map(cvcl2name).fillna(K.line.astype(str))
    K["k"] = K.drug.map(norm)
    out = {}
    for k, gg in K.groupby("k", observed=True):
        by_line = {}
        for cn, gl in gg.groupby("cname", observed=True):
            by_line[norm(cn)] = DELTA[gl.i.to_numpy()].mean(0)
        if len(by_line) < 2:
            continue
        names = list(by_line)
        A = np.stack([by_line[c] for c in names])
        tot = A.sum(0)
        for j, c in enumerate(names):
            out[(c, k)] = A[j] - (tot - A[j]) / (len(names) - 1)
    # remove each line's general response before calling the rest a
    # cell-drug relation -- see celldrug.remove_line_effect_profiles
    return remove_line_effect_profiles(out), genes


def identity_check(a, ga, b, gb, label, min_cpd=8):
    """Match cell lines by DATA, not by name.

    Names are a weak basis for saying two datasets measured the same line:
    they disagree in formatting (MCF-7 / MCF7 / MCF7.311), and Ben-David et al.
    showed the same nominal line genuinely diverges between laboratories, so a
    name match is a hypothesis rather than a fact.

    Baseline-expression matching would be the natural check but is impossible
    here: LINCS Level 4 is z-scored within plate, which removes exactly the
    cell-line baseline identity we would want to match on. What survives that
    normalisation is the RESPONSE, so we match on the response fingerprint --
    each line's residual profile concatenated across the compounds both datasets
    share. If two datasets really measured the same line, that line should be
    its own best match across all candidates.

    Returns the similarity matrix and a per-pair table with the rank a
    name-matched partner achieves, so name matches can be kept, rejected, or
    replaced by the reciprocal best hit.
    """
    ia = {g: i for i, g in enumerate(ga)}
    shared_g = [g for g in gb if g in ia]
    if len(shared_g) < MIN_GENES:
        return None, []
    ai = np.array([ia[g] for g in shared_g])
    ib = {g: i for i, g in enumerate(gb)}
    bi = np.array([ib[g] for g in shared_g])
    la = sorted({k[0] for k in a}); lb = sorted({k[0] for k in b})
    cpd_a = {k[1] for k in a}; cpd_b = {k[1] for k in b}
    shared_c = sorted(cpd_a & cpd_b)
    if len(shared_c) < min_cpd:
        print(f"  {label}: only {len(shared_c)} shared compounds — "
              f"identity check skipped")
        return None, []

    def fingerprint(store, lines, idx):
        F = {}
        for ln in lines:
            vecs, used = [], []
            for c in shared_c:
                v = store.get((ln, c))
                if v is not None:
                    vecs.append(v[idx]); used.append(c)
            if len(used) >= min_cpd:
                F[ln] = (dict(zip(used, vecs)))
        return F

    FA, FB = fingerprint(a, la, ai), fingerprint(b, lb, bi)
    rows = []
    M = pd.DataFrame(np.nan, index=sorted(FA), columns=sorted(FB))
    for x in FA:
        for y in FB:
            common = set(FA[x]) & set(FB[y])
            if len(common) < min_cpd:
                continue
            u = np.concatenate([FA[x][c] for c in sorted(common)])
            v = np.concatenate([FB[y][c] for c in sorted(common)])
            if np.std(u) > 0 and np.std(v) > 0:
                M.loc[x, y] = float(stats.pearsonr(u, v).statistic)
    for x in M.index:
        r = M.loc[x].dropna()
        if not len(r):
            continue
        order = r.sort_values(ascending=False)
        best = order.index[0]
        rank = int(list(order.index).index(x) + 1) if x in order.index else -1
        rows.append({"comparison": label, "line": x, "name_match_present":
                     x in order.index,
                     "r_name_match": float(order.get(x, np.nan)),
                     "best_match": best, "r_best": float(order.iloc[0]),
                     "rank_of_name_match": rank, "n_candidates": len(order),
                     "reciprocal_best": bool(best == x)})
    D = pd.DataFrame(rows)
    if len(D):
        nm = D[D.name_match_present]
        print(f"  {label} IDENTITY: {len(nm)} name-matched lines; "
              f"{int(nm.reciprocal_best.sum())} are their own best match "
              f"({nm.reciprocal_best.mean():.0%}); median rank of the name "
              f"match {nm.rank_of_name_match.median():.0f} of "
              f"{nm.n_candidates.median():.0f}")
    return M, rows


def load_sciplex(path):
    """sciPlex3 (Srivatsan et al., Trapnell lab) as a THIRD laboratory.

    Tahoe and LINCS are two institutions; a third is what turns a comparison
    into a trend. sciPlex3 shares A549, K562 and MCF7 with LINCS and uses a
    different technology again (nuclear hashing rather than L1000), so
    agreement across it is not attributable to a shared platform.
    """
    import anndata as ad
    a = ad.read_h5ad(path)
    obs = a.obs.copy()
    obs["k"] = obs.perturbation.astype(str).map(norm)
    obs["cl"] = obs.cell_line.astype(str).map(norm)
    obs = obs[(obs.cl != "nan") & (obs.k != "nan") & (obs.k != "")]
    X = np.asarray(a.X.todense()) if hasattr(a.X, "todense") else np.asarray(a.X)
    X = X[obs.index.map(lambda i: a.obs.index.get_loc(i)).to_numpy()]
    tot = X.sum(1, keepdims=True)
    X = np.log1p(X / np.where(tot > 0, tot, 1) * 1e6).astype(np.float32)
    genes = np.array([str(g) for g in a.var_names])
    # Subtract each line's own vehicle. Without it the "residual" is baseline
    # expression: measured against the file's own control wells, the previous
    # residuals correlated with the same line's vehicle at median r = 0.898.
    veh = {}
    for cl, gl in obs[obs.k.isin({"control", "vehicle", "dmso"})].groupby(
            "cl", observed=True):
        ii = [obs.index.get_loc(i) for i in gl.index]
        veh[cl] = X[ii].mean(0)
    if not veh:
        print("  sciPlex3: no vehicle wells found; arm skipped", flush=True)
        return {}, genes
    out = {}
    for k, gg in obs.groupby("k", observed=True):
        if k in {"control", "vehicle", "dmso"}:
            continue
        by_line = {}
        for cl, gl in gg.groupby("cl", observed=True):
            if cl not in veh:
                continue
            ii = [obs.index.get_loc(i) for i in gl.index]
            by_line[cl] = X[ii].mean(0) - veh[cl]
        if len(by_line) < 2:
            continue
        names = list(by_line)
        A = np.stack([by_line[c] for c in names])
        tot_ = A.sum(0)
        for j, c in enumerate(names):
            out[(c, k)] = A[j] - (tot_ - A[j]) / (len(names) - 1)
    print(f"  sciPlex3: {len({k[0] for k in out})} lines, "
          f"{len({k[1] for k in out})} compounds", flush=True)
    # remove each line's general response before calling the rest a
    # cell-drug relation -- see celldrug.remove_line_effect_profiles
    return remove_line_effect_profiles(out), genes


def compare(a, ga, b, gb, label, allowed=None):
    """Correlate residual profiles on the shared genes, per (line, compound)."""
    ia = {g: i for i, g in enumerate(ga)}
    shared = [g for g in gb if g in ia]
    if len(shared) < MIN_GENES:
        print(f"  {label}: only {len(shared)} shared genes — skipped")
        return []
    ai = np.array([ia[g] for g in shared])
    ib = {g: i for i, g in enumerate(gb)}
    bi = np.array([ib[g] for g in shared])
    keys = set(a) & set(b)
    if allowed is not None:
        keys = {k for k in keys if k[0] in allowed}
    rows = []
    for (line, cpd) in keys:
        x = a[(line, cpd)][ai]
        y = b[(line, cpd)][bi]
        if np.std(x) == 0 or np.std(y) == 0:
            continue
        r = stats.pearsonr(x, y)
        rows.append({"comparison": label, "line": line, "compound": cpd,
                     "n_genes": len(shared), "r": float(r.statistic),
                     "p": float(r.pvalue)})
    print(f"  {label}: {len(rows)} (line, compound) pairs over "
          f"{len(shared)} shared genes", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    for gz, gc in ((LIN / "level4.gctx.gz", LIN / "level4.gctx"),
                   (LIN / "p1_level4.gctx.gz", LIN / "p1_level4.gctx")):
        if not gc.exists() and gz.exists():
            print(f"decompressing {gz.name} ...", flush=True)
            with gzip.open(gz, "rb") as f, open(gc, "wb") as o:
                shutil.copyfileobj(f, o, length=1 << 24)

    print("loading LINCS phase 1 ...", flush=True)
    p1, g1 = load_lincs(LIN / "p1_level4.gctx", LIN / "p1_inst_info.txt.gz",
                        LIN / "p1_gene_info.txt.gz", "phase1")
    print("loading LINCS phase 2 ...", flush=True)
    p2, g2 = load_lincs(LIN / "level4.gctx", LIN / "inst_info.txt.gz",
                        LIN / "GSE70138_Broad_LINCS_gene_info.txt.gz", "phase2")

    ident = []
    M_wl, r_wl = identity_check(p1, g1, p2, g2,
                                "LINCS p1 vs p2 (within-lab, Broad)")
    ident += r_wl
    ok_wl = {r["line"] for r in r_wl if r["reciprocal_best"]} or None
    rows = compare(p1, g1, p2, g2, "LINCS p1 vs p2 (within-lab, Broad)")
    rows += compare(p1, g1, p2, g2,
                    "LINCS p1 vs p2 (within-lab, identity-validated)",
                    allowed=ok_wl)

    print("loading Tahoe ...", flush=True)
    try:
        th, gt = load_tahoe(args.pb_dir)
        for store, gg, tag in ((p1, g1, "p1"), (p2, g2, "p2")):
            Mx, rx = identity_check(th, gt, store, gg,
                                    f"Tahoe vs LINCS {tag} (CROSS-LAB)")
            ident += rx
            okx = {r["line"] for r in rx if r["reciprocal_best"]} or None
            rows += compare(th, gt, store, gg,
                            f"Tahoe vs LINCS {tag} (CROSS-LAB)")
            rows += compare(th, gt, store, gg,
                            f"Tahoe vs LINCS {tag} (CROSS-LAB, "
                            f"identity-validated)", allowed=okx)
    except Exception as e:
        print(f"  Tahoe arm unavailable: {e}")

    sp_path = ROOT / "data/external/scperturb/sciplex3_pseudobulk.h5ad"
    if sp_path.exists():
        print("loading sciPlex3 (third laboratory) ...", flush=True)
        try:
            sp, gsp = load_sciplex(sp_path)
            for store, gg, tag in ((p1, g1, "p1"), (p2, g2, "p2")):
                rows += compare(sp, gsp, store, gg,
                                f"sciPlex3 vs LINCS {tag} (CROSS-LAB)")
        except Exception as e:
            print(f"  sciPlex3 arm unavailable: {e}")

    R = pd.DataFrame(rows)
    if not len(R):
        print("no comparisons produced"); return
    R.to_csv(TAB / "cross_lab_transcription.csv", index=False)
    if ident:
        I = pd.DataFrame(ident)
        I.to_csv(TAB / "cross_lab_identity.csv", index=False)
        print("\n=== cell-line identity by response fingerprint ===")
        for cmpn, g in I.groupby("comparison"):
            nm = g[g.name_match_present]
            if not len(nm):
                print(f"  {cmpn}: no name-matched lines in common")
                continue
            print(f"  {cmpn}: {int(nm.reciprocal_best.sum())}/{len(nm)} "
                  f"name matches are reciprocal best hits; median r "
                  f"name-match {nm.r_name_match.median():.3f} vs best "
                  f"{nm.r_best.median():.3f}")
        print("  A name match that is NOT its own best hit means the two "
              "datasets'\n  nominally identical lines behave less like each "
              "other than like some\n  other line — which is what Ben-David "
              "predicts and what name-based\n  matching silently assumes away.")
    S = R.groupby("comparison").r.agg(["median", "mean", "size",
                                       lambda s: float((s > 0).mean())])
    S.columns = ["median_r", "mean_r", "n_pairs", "frac_positive"]
    print("\n=== residual-profile correlation per (line, compound) ===")
    print(S.round(3).to_string())

    # The headline fraction must be built the same way as the viability arm:
    # non-identity-validated rows only, and excluding sciPlex3, whose three
    # contexts make its leave-one-context-out residual a different quantity
    # (see RESULTS.md 23). Pooling those rows in silently moved this number.
    # Use ONE cross-lab comparison, not the pooled pair: 228 of the 489 pooled
    # rows were the same (line, compound) counted twice, once against phase 1
    # and once against phase 2. And measure the ceiling on the SAME (line,
    # compound) pairs as the numerator, which the previous ratio did not --
    # it divided 489 cross-lab pairs over 3-6 lines by 5,803 within-lab pairs
    # over 16 lines.
    wl_all = R[R.comparison.str.contains("within-lab")
               & ~R.comparison.str.contains("validated")]
    xl = R[(R.comparison == "Tahoe vs LINCS p1 (CROSS-LAB)")]
    keys = set(zip(xl.line, xl.compound))
    wl_m = wl_all[[(l, c) in keys for l, c in zip(wl_all.line, wl_all.compound)]]
    within = float(wl_m.r.median()) if len(wl_m) >= 20 else float(wl_all.r.median())
    cross = float(xl.r.median()) if len(xl) else np.nan
    print(f"\nmatched ceiling: {len(wl_m)} within-lab pairs sharing a "
          f"(line, compound) with the {len(xl)} cross-lab pairs", flush=True)
    if np.isfinite(within) and np.isfinite(cross) and within > 0:
        print(f"\nwithin-lab (Broad, p1 vs p2): {within:.3f}")
        print(f"cross-lab (Tahoe vs LINCS):   {cross:.3f}")
        print(f"REPRODUCIBLE FRACTION (transcription) = {cross/within:.1%}")
        print("  Compare the viability arm's 56%. Agreement between the two "
              "readouts would\n  make this a property of laboratories rather "
              "than of a particular assay.")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    order = list(S.index)
    cols = [AQUA if "within" in o else ORANGE for o in order]
    data = [R[R.comparison == o].r.dropna() for o in order]
    bp = ax[0].boxplot(data, showfliers=False, patch_artist=True,
                       medianprops=dict(color="black", lw=1.6))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.75)
    ax[0].set_xticks(range(1, len(order) + 1),
                     [o.split("(")[0].strip()[:18] for o in order],
                     fontsize=7, rotation=15, ha="right")
    ax[0].axhline(0, color="#444444", lw=0.9)
    ax[0].set_ylabel("residual-profile r")
    for i, d in enumerate(data):
        if len(d):
            ax[0].text(i + 1, np.median(d) + 0.01, f"{np.median(d):.2f}",
                       ha="center", fontsize=9, fontweight="bold")
    ax[0].set_title("A  Transcriptional residual reproducibility", loc="left",
                    fontweight="bold", fontsize=10)
    for o, c in zip(order, cols):
        d = R[R.comparison == o].r.dropna()
        if len(d):
            ax[1].hist(d, bins=40, histtype="step", lw=2, color=c, density=True,
                       label=o.split("(")[0].strip()[:22])
    ax[1].axvline(0, color="#444444", lw=0.9)
    ax[1].set_xlabel("residual-profile r"); ax[1].set_ylabel("density")
    ax[1].legend(frameon=False, fontsize=7)
    ax[1].set_title("B  Distributions", loc="left", fontweight="bold",
                    fontsize=10)
    fig.suptitle("Does the line-specific TRANSCRIPTIONAL response reproduce "
                 "across laboratories?", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "cross_lab_transcription", FIG,
                    source_data={"per_pair": R, "summary": S.reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
