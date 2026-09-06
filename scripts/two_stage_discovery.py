#!/usr/bin/env python3
"""Detect the interaction first, then ask what explains it.

RESULTS.md sec.37 found that a genome-wide ridge over 11,809 mutated genes has no
predictive power for drug response (CV R^2 = -0.019) while a single known marker
reaches +0.080 -- the signal is sparse and averaging destroys it. Sec.39 then
built a detector that finds WHICH (line, drug) pairs carry reproducible
interaction, without prior knowledge.

Putting the two together gives a method that neither half supplies alone:

  STAGE 1  Detect. Per (line, compound), the dose-response residual is measured
           on replicate detection plate X1 and again on X2/X3. Their agreement
           says whether that pair's interaction is real. No marker, pathway or
           genotype is used.
  STAGE 2  Associate, but ONLY where stage 1 found something. A marker scan run
           over every pair spends its multiple-testing budget on pairs that
           carry no interaction at all; restricted to reproducible pairs, the
           same scan tests far fewer hypotheses and each has a real effect to
           find.

This is the standard two-stage logic of filtering on an independent statistic
before testing, and the independence here is structural: stage 1 uses only
cross-replicate agreement and never looks at any molecular feature, so it cannot
bias the stage-2 association.

Three things are asked of it:

  BENCHMARK   Does filtering recover the established pharmacogenomic
              relationships of sec.37 -- MDM2/TP53, BRAF, MEK/RAS, EGFR -- better
              than the same scan unfiltered? If not, the method is not worth
              having.
  DISCOVERY   What associations does it find that are NOT already-known
              pharmacology?
  VALIDATION  Do those hold in GDSC, a different laboratory with different
              chemistry? PRISM and GDSC share cell lines and compounds, so a
              real drug x genotype relationship should appear in both.

Outputs: results/tables/two_stage_discovery.csv
         results/tables/two_stage_benchmark.csv
         figure bundle results/figures/00_manuscript/two_stage_discovery/
"""
import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.sparse_interaction import _derangement
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PR = ROOT / "data" / "external" / "prism"
GD = ROOT / "data" / "external" / "gdsc"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"

# the established relationships of sec.37, used here only to BENCHMARK the
# method -- they are never supplied to it
KNOWN = {"TP53": ["nutlin-3", "idasanutlin", "amg-232"],
         "BRAF": ["vemurafenib", "dabrafenib", "plx-4720"],
         "KRAS": ["trametinib", "cobimetinib", "selumetinib"],
         "NRAS": ["trametinib", "cobimetinib", "selumetinib"],
         "EGFR": ["erlotinib", "gefitinib", "osimertinib"]}


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=300)
    ap.add_argument("--min-lines", type=int, default=25)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    lfc = pd.read_csv(PR / "secondary-screen-logfold-change.csv", index_col=0)
    lfc = lfc[[not i.endswith("_FAILED_STR") for i in lfc.index]]
    lfc.index = [re.search(r"(ACH-\d+)", i).group(1) for i in lfc.index]
    lfc = lfc.groupby(level=0).mean()
    ti = pd.read_csv(PR / "secondary-screen-replicate-treatment-info.csv",
                     low_memory=False)
    ti = ti[ti.column_name.isin(lfc.columns) & ti.name.notna()].copy()
    ti["rep"] = ti.detection_plate.astype(str).str.extract(r"_(X\d)")[0]
    ti = ti[ti.rep.notna()]
    ti["dose_s"] = ti.dose.round(4).astype(str)
    L = lfc.to_numpy(np.float32)
    lines = np.array(lfc.index)
    colpos = {c: i for i, c in enumerate(lfc.columns)}
    print(f"PRISM {L.shape[0]} lines x {ti.name.nunique()} compounds",
          flush=True)

    # ---------- STAGE 1: detect, using replicate plates only ----------
    # dose-response residual per (line, compound, replicate half)
    # Compounds are screened at different numbers of doses, so the residual
    # vectors have different lengths and cannot go into one array. Fix the
    # length: use each compound's own first N_DOSE doses, sorted numerically,
    # and skip compounds with fewer. Sorting the dose STRINGS would order
    # "10" before "2", which would scramble the dose axis the cosine is
    # computed along.
    N_DOSE = 6
    A, B, meta = [], [], []
    for cpd, gc in ti.groupby("name", observed=True):
        doses = [d for _, d in sorted(
            {(float(x), x) for x in gc.dose_s.unique()})][:N_DOSE]
        if len(doses) < N_DOSE:
            continue
        halves = {}
        for tag, reps in (("A", ["X1"]), ("B", ["X2", "X3"])):
            mats = []
            for d in doses:
                sub = gc[(gc.dose_s == d) & gc.rep.isin(reps)]
                if not len(sub):
                    mats.append(np.full(L.shape[0], np.nan)); continue
                idx = [colpos[c] for c in sub.column_name]
                with np.errstate(invalid="ignore"):
                    mats.append(np.nanmean(L[:, idx], axis=1))
            halves[tag] = np.stack(mats, axis=1)          # lines x doses
        Ha, Hb = halves["A"], halves["B"]
        # residual: remove the compound's leave-one-line-out mean at each dose
        for H in (Ha, Hb):
            with np.errstate(invalid="ignore"):
                tot = np.nansum(H, axis=0, keepdims=True)
                cnt = np.sum(np.isfinite(H), axis=0, keepdims=True)
                loo = (tot - np.nan_to_num(H)) / np.maximum(cnt - 1, 1)
            H -= loo
        ok = np.isfinite(Ha).all(1) & np.isfinite(Hb).all(1)
        if ok.sum() < args.min_lines:
            continue
        for i in np.where(ok)[0]:
            A.append(Ha[i]); B.append(Hb[i]); meta.append((lines[i], cpd))
    A, B = np.stack(A), np.stack(B)
    M = pd.DataFrame(meta, columns=["line", "compound"])
    print(f"stage 1: {len(M)} (line, compound) pairs with a dose-response "
          f"vector on both replicate halves", flush=True)

    # remove each line's general sensitivity, leave-one-compound-out
    for tag, Z in (("A", A), ("B", B)):
        for ln, g in M.groupby("line", observed=True):
            ii = g.index.to_numpy()
            if len(ii) < 5:
                continue
            tot = Z[ii].sum(0)
            Z[ii] -= (tot - Z[ii]) / (len(ii) - 1)

    na, nb = np.linalg.norm(A, axis=1), np.linalg.norm(B, axis=1)
    ok = (na > 0) & (nb > 0)
    cos = np.full(len(A), np.nan)
    cos[ok] = (A[ok] * B[ok]).sum(1) / (na[ok] * nb[ok])
    rng = np.random.default_rng(0)
    blk = M.line.to_numpy()
    null = []
    for _ in range(args.n_perm):
        idx = _derangement(int(ok.sum()), rng, blk[ok])
        null.append((A[ok] * B[ok][idx]).sum(1) / (na[ok] * nb[ok][idx]))
    flat = np.sort(np.concatenate(null))
    p = np.full(len(A), np.nan)
    p[ok] = 1.0 - np.searchsorted(flat, cos[ok]) / len(flat)
    p[ok] = np.clip(p[ok], 1.0 / (len(flat) + 1), 1.0)
    order = np.argsort(np.where(np.isfinite(p), p, 2))
    q = np.full(len(p), np.nan)
    prev = 1.0
    n_ok = int(ok.sum())
    for r_, i in enumerate(order[:n_ok][::-1]):
        prev = min(prev, p[i] * n_ok / (n_ok - r_))
        q[i] = prev
    M["cos"], M["p"], M["q"] = cos, p, q
    flagged = M.q < 0.05
    print(f"   {int(flagged.sum())} of {n_ok} pairs reproduce at FDR<0.05 "
          f"({flagged.sum()/n_ok:.1%})", flush=True)

    # ---------- molecular features ----------
    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    dep = smp.cross_references.astype(str).str.extract(
        r"DepMap;\s*(ACH-\d+)")[0]
    cv2dep = dict(zip(smp.Cellosaurus_ID, dep))
    mut = pd.read_csv(PR / "mutation_long.tsv.gz", sep="\t", low_memory=False)
    mut["depmap"] = mut.RRID.map(cv2dep)
    ns = mut[mut.depmap.notna() &
             mut.Variant_Classification.astype(str).str.contains(
                 "Missense|Nonsense|Frame_Shift|Splice|In_Frame", na=False)]
    called = set(ns.depmap)
    gc_ = ns.groupby("Gene_symbol").depmap.nunique()
    genes = sorted(gc_[gc_ >= 25].index)
    gsets = {g: set(d.depmap) for g, d in
             ns[ns.Gene_symbol.isin(genes)].groupby("Gene_symbol")}
    print(f"   {len(genes):,} genes mutated in >=25 sequenced lines", flush=True)

    # ---------- STAGE 2: associate, filtered vs unfiltered ----------
    def scan(sub, tag):
        rows = []
        for cpd, g in sub.groupby("compound", observed=True):
            ls = [l for l in g.line if l in called]
            if len(ls) < args.min_lines:
                continue
            y = g.set_index("line").cos.reindex(ls).to_numpy()
            for gene in genes:
                m = np.array([l in gsets[gene] for l in ls])
                if m.sum() < 5 or (~m).sum() < 10:
                    continue
                u = stats.mannwhitneyu(y[m], y[~m], alternative="two-sided")
                rows.append({"compound": cpd, "gene": gene, "n_mut": int(m.sum()),
                             "n_wt": int((~m).sum()),
                             "delta": float(np.median(y[m]) - np.median(y[~m])),
                             "p": float(u.pvalue), "scan": tag})
        R = pd.DataFrame(rows)
        if len(R):
            R = R.sort_values("p")
            R["q"] = np.minimum(R.p * len(R) / np.arange(1, len(R) + 1), 1)
            R["q"] = R.q[::-1].cummin()[::-1]
        return R

    print("\nstage 2: marker scan, filtered vs unfiltered", flush=True)
    F = scan(M[flagged], "filtered")
    U = scan(M, "unfiltered")
    print(f"   filtered   {len(F):,} tests, "
          f"{int((F.q < 0.05).sum()) if len(F) else 0} at FDR<0.05")
    print(f"   unfiltered {len(U):,} tests, "
          f"{int((U.q < 0.05).sum()) if len(U) else 0} at FDR<0.05")

    # ---------- BENCHMARK on the known relationships ----------
    def rank_of_known(R):
        out = []
        if not len(R):
            return pd.DataFrame()
        R = R.reset_index(drop=True)
        for gene, drugs in KNOWN.items():
            for d in drugs:
                hit = R[(R.gene == gene)
                        & (R.compound.map(norm) == norm(d))]
                if len(hit):
                    out.append({"gene": gene, "compound": d,
                                "p": float(hit.p.iloc[0]),
                                "q": float(hit.q.iloc[0]),
                                "rank": int(hit.index[0]) + 1,
                                "n_tests": len(R)})
        return pd.DataFrame(out)

    BF, BU = rank_of_known(F), rank_of_known(U)
    print("\nBENCHMARK — where the known relationships rank in each scan")
    for tag, B_ in (("filtered", BF), ("unfiltered", BU)):
        if len(B_):
            print(f"   {tag:11s} {len(B_)} testable; median rank "
                  f"{B_['rank'].median():.0f} of {int(B_.n_tests.iloc[0]):,}; "
                  f"{int((B_.q < 0.05).sum())} at FDR<0.05")
        else:
            print(f"   {tag:11s} none testable")
    if len(BF) and len(BU):
        j = BF.merge(BU, on=["gene", "compound"], suffixes=("_f", "_u"))
        if len(j):
            better = int((j.rank_f / j.n_tests_f < j.rank_u / j.n_tests_u).sum())
            print(f"   the known relationship ranks relatively higher in the "
                  f"filtered scan for {better}/{len(j)}")

    # ---------- DISCOVERY: hits that are not known pharmacology -------
    known_pairs = {(g, norm(d)) for g, ds in KNOWN.items() for d in ds}
    if len(F):
        top = F[F.q < 0.05].copy()
        top["known"] = [(r.gene, norm(r.compound)) in known_pairs
                        for r in top.itertuples()]
        nov = top[~top.known]
        print(f"\nDISCOVERY — {len(top)} associations at FDR<0.05, "
              f"{len(nov)} not in the benchmark set")
        if len(nov):
            print(nov.head(12)[["compound", "gene", "n_mut", "delta",
                                "q"]].round(4).to_string(index=False))
        F.to_csv(TAB / "two_stage_discovery.csv", index=False)
    else:
        nov = pd.DataFrame()

    # ---------- VALIDATION in GDSC ----------
    print("\nVALIDATION — do the new associations hold in GDSC?", flush=True)
    val = pd.DataFrame()
    try:
        cos_map = {}
        for d, xr in zip(dep, smp.cross_references.astype(str)):
            if isinstance(d, str):
                for cid in re.findall(r"Cosmic(?:-CLP)?;\s*(\d+)", xr):
                    cos_map.setdefault(int(cid), d)
        gd = []
        for f in ("GDSC1_fitted_dose_response_27Oct23.csv",
                  "GDSC2_fitted_dose_response_27Oct23.csv"):
            if (GD / f).exists():
                gd.append(pd.read_csv(GD / f, low_memory=False))
        Gd = pd.concat(gd, ignore_index=True)
        Gd["dep"] = Gd.COSMIC_ID.map(cos_map)
        Gd = Gd[Gd.dep.notna()]
        Gd["k"] = Gd.DRUG_NAME.map(norm)
        rows = []
        for r_ in (nov.head(40).itertuples() if len(nov) else []):
            sub = Gd[Gd.k == norm(r_.compound)]
            if len(sub) < 40:
                continue
            v = sub.groupby("dep").Z_SCORE.mean()
            m = np.array([d in gsets[r_.gene] for d in v.index])
            if m.sum() < 5 or (~m).sum() < 15:
                continue
            u = stats.mannwhitneyu(v[m], v[~m], alternative="two-sided")
            rows.append({"compound": r_.compound, "gene": r_.gene,
                         "prism_delta": r_.delta, "prism_q": r_.q,
                         "gdsc_delta": float(v[m].median() - v[~m].median()),
                         "gdsc_p": float(u.pvalue), "n_gdsc_mut": int(m.sum())})
        val = pd.DataFrame(rows)
        if len(val):
            val["q"] = np.minimum(val.gdsc_p * len(val), 1)
            n_v = int((val.q < 0.05).sum())
            print(f"   {len(val)} of the new associations are testable in GDSC; "
                  f"{n_v} replicate at Bonferroni q<0.05")
            print(val.sort_values("gdsc_p").head(10).round(4).to_string(
                index=False))
            val.to_csv(TAB / "two_stage_validation.csv", index=False)
        else:
            print("   none testable in GDSC")
    except Exception as e:
        print(f"   GDSC validation unavailable: {type(e).__name__}: {e}")

    bench = pd.concat([BF.assign(scan="filtered"),
                       BU.assign(scan="unfiltered")], ignore_index=True) \
        if len(BF) or len(BU) else pd.DataFrame()
    if len(bench):
        bench.to_csv(TAB / "two_stage_benchmark.csv", index=False)

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3), constrained_layout=True)
    ax[0].hist(M.cos.dropna(), bins=60, color=GREY, alpha=0.85)
    if flagged.any():
        ax[0].hist(M[flagged].cos, bins=30, color=ORANGE, alpha=0.9,
                   label=f"{int(flagged.sum())} at FDR<0.05")
        ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_xlabel("cross-replicate agreement of the dose-response residual")
    ax[0].set_ylabel("(line, compound) pairs")
    ax[0].set_title("a  Stage 1: detect, without markers", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(bench):
        d0 = [bench[bench.scan == s]["rank"] / bench[bench.scan == s].n_tests
              for s in ("filtered", "unfiltered")
              if (bench.scan == s).any()]
        if len(d0) == 2:
            bp = ax[1].boxplot(d0, showfliers=False, patch_artist=True,
                               medianprops=dict(color="black", lw=1.4))
            for pch, c in zip(bp["boxes"], [ORANGE, GREY]):
                pch.set_facecolor(c); pch.set_alpha(0.85)
            ax[1].set_xticks([1, 2], ["filtered", "unfiltered"], fontsize=8)
            ax[1].set_yscale("log")
            ax[1].set_ylabel("relative rank of the known relationship")
    ax[1].set_title("b  Benchmark on known pharmacology", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(val):
        ax[2].scatter(val.prism_delta, val.gdsc_delta, s=40,
                      color=[ORANGE if q < 0.05 else GREY for q in val.q],
                      edgecolors="none")
        ax[2].axhline(0, color="#888", lw=0.8); ax[2].axvline(0, color="#888",
                                                              lw=0.8)
        r_ = stats.spearmanr(val.prism_delta, val.gdsc_delta)
        ax[2].set_xlabel("effect in PRISM (discovery)")
        ax[2].set_ylabel("effect in GDSC (validation)")
        ax[2].text(0.04, 0.94, f"ρ = {r_.statistic:+.2f}\n"
                   f"{int((val.q < 0.05).sum())}/{len(val)} replicate",
                   transform=ax[2].transAxes, va="top", fontsize=8)
    ax[2].set_title("c  Independent laboratory", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Two-stage discovery: detect the interaction first, then ask "
                 "what explains it", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    d = save_figure(fig, "two_stage_discovery", FIG,
                    source_data={"pairs": M, "filtered": F,
                                 "benchmark": bench, "validation": val},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
