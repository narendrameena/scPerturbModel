#!/usr/bin/env python3
"""Independent validation of the allele candidates in GDSC.

The four "novel" allele x compound associations (RESULTS.md §15) passed
split-sample replication at 94-99%, but that test is circular: candidates were
selected using every cell line, so re-splitting those same lines measures effect
size, not generalisation. This script does the non-circular version -- a
different laboratory, a different assay (fitted IC50 rather than pooled
barcode viability), and cell lines whose response was never used to pick the
candidates.

**The pipeline has to prove it works before its negatives mean anything.** A
validation that finds nothing is uninformative if it would also have found
nothing for a true effect. Every run therefore tests known biomarkers first --
BRAF V600E against dabrafenib/trametinib/selumetinib, PIK3CA hotspots against
alpelisib -- and the candidate results are only interpretable if those pass.

**Compound coverage is the limiting factor and is reported, not hidden.** Of the
four candidates only lapatinib exists in GDSC; poziotinib, BMS-690514 and
mozavaptan do not. So each candidate is tested two ways:

  exact      the same compound, where GDSC has it
  class      every GDSC compound sharing the candidate's mechanism. Three of the
             four candidates are EGFR/HER inhibitors, and the claim being made
             is that the allele modulates response to that class -- which is
             testable against five EGFR inhibitors even when the exact molecule
             is missing, and is a fairer test of the biology than one compound.

Direction matters and is checked, not just significance. GDSC Z_SCORE is
z-scored log IC50 within drug, so lower = more sensitive; PRISM effects are
residual log-fold-change, where negative = more killing = more sensitive. A real
association must reproduce with the SAME sign.

Outputs: results/tables/gdsc_allele_validation.csv
         figure bundle results/figures/13_prism/gdsc_validation/
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

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
GD = ROOT / "data" / "external" / "gdsc"
PR = ROOT / "data" / "external" / "prism"
FIG = ROOT / "results" / "figures" / "13_prism"
TAB = ROOT / "results" / "tables"
MIN_MUT = 5
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"

# candidate -> (allele, PRISM effect sign, mechanism keywords for the class test)
CANDIDATES = [
    ("mozavaptan", "ZMYM5 p.K463fs", -1, ["AVPR1", "AVPR2", "vasopressin"]),
    ("BMS-690514", "MUC16 p.T12415A", -1, ["EGFR", "ERBB2", "ERBB3", "HER2"]),
    ("poziotinib", "HMCN1 p.K2374fs", +1, ["EGFR", "ERBB2", "ERBB3", "HER2"]),
    ("lapatinib", "MUC16 p.T12415A", -1, ["EGFR", "ERBB2", "ERBB3", "HER2"]),
]
# known biomarkers used as positive controls: (allele, drugs, expected sign)
CONTROLS = [
    ("BRAF p.V600E", ["dabrafenib", "trametinib", "selumetinib", "PLX-4720",
                      "SB590885"], -1),
    ("PIK3CA p.E545K", ["alpelisib", "pictilisib", "taselisib"], -1),
    ("PIK3CA p.H1047R", ["alpelisib", "pictilisib", "taselisib"], -1),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metric", default="Z_SCORE",
                    choices=["Z_SCORE", "LN_IC50", "AUC"])
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    g1 = pd.read_csv(GD / "GDSC1_fitted_dose_response_27Oct23.csv",
                     low_memory=False)
    g2 = pd.read_csv(GD / "GDSC2_fitted_dose_response_27Oct23.csv",
                     low_memory=False)
    G = pd.concat([g1.assign(SET="GDSC1"), g2.assign(SET="GDSC2")],
                  ignore_index=True)
    G["drug_l"] = G.DRUG_NAME.astype(str).str.lower().str.strip()
    print(f"GDSC: {len(G)} measurements, {G.DRUG_NAME.nunique()} drugs, "
          f"{G.COSMIC_ID.nunique()} cell lines", flush=True)

    # COSMIC id -> Cellosaurus RRID, from the Cellosaurus cross-reference field
    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    cos2rrid = {}
    for rrid, xr in zip(smp.Cellosaurus_ID, smp.cross_references.astype(str)):
        for cid in re.findall(r"Cosmic(?:-CLP)?;\s*(\d+)", xr):
            cos2rrid.setdefault(int(cid), rrid)
    G["RRID"] = G.COSMIC_ID.map(cos2rrid)
    mapped = G.RRID.notna()
    print(f"mapped {G.loc[mapped, 'COSMIC_ID'].nunique()} of "
          f"{G.COSMIC_ID.nunique()} GDSC cell lines to Cellosaurus IDs "
          f"({mapped.mean():.0%} of measurements)", flush=True)
    G = G[mapped]

    lineage = (G.drop_duplicates("RRID").set_index("RRID").TCGA_DESC
               .astype(str))
    mut = pd.read_csv(PR / "mutation_long.tsv.gz", sep="\t", low_memory=False)
    tmb = mut.groupby("RRID").size()          # total mutational burden proxy
    dmg = mut[mut.Variant_Classification.astype(str).str.contains(
        "Missense|Nonsense|Frame_Shift|Splice|In_Frame", na=False)].copy()
    dmg["pc"] = dmg.Protein_Change.astype(str).str.strip()
    dmg["vid"] = dmg.Gene_symbol.astype(str) + " " + dmg.pc
    vsets = {v: set(d.RRID) for v, d in dmg.groupby("vid")}
    dmg_all = mut.copy()          # incl. synonymous, for the germline check
    gsets = {g: set(d.RRID) for g, d in dmg.groupby("Gene_symbol")}

    def test(allele, drugs, expect_sign, label, kind):
        rows = []
        for dn in drugs:
            sub = G[G.drug_l == dn.lower()]
            if not len(sub):
                continue
            # one value per cell line (GDSC1/GDSC2 and repeats averaged)
            v = sub.groupby("RRID")[args.metric].mean()
            carriers = vsets.get(allele, set())
            m = v.index.isin(carriers)
            if m.sum() < MIN_MUT or (~m).sum() < MIN_MUT:
                rows.append({"label": label, "kind": kind, "allele": allele,
                             "drug": dn, "n_mut": int(m.sum()),
                             "n_wt": int((~m).sum()), "delta": np.nan,
                             "p": np.nan, "expect_sign": expect_sign,
                             "note": "too few carriers in GDSC"})
                continue
            u, p = stats.mannwhitneyu(v[m], v[~m], alternative="two-sided")
            d = float(np.median(v[m]) - np.median(v[~m]))
            # Two confounds that would manufacture exactly this result:
            #  lineage -- alleles are unevenly distributed across tissues, and
            #             EGFR-inhibitor sensitivity is strongly lineage-driven
            #  TMB     -- MUC16 is the classic tumour-mutational-burden proxy
            #             (2nd most mutated gene of 18,739), so "MUC16 mutant"
            #             is close to "highly mutated line"
            # Both are removed by residualising the response on them before
            # testing, rank-based so no distributional assumption is added.
            lin = lineage.reindex(v.index).fillna("?")
            zz = v.groupby(lin).rank(pct=True)          # within-lineage rank
            p_lin = (stats.mannwhitneyu(zz[m], zz[~m])[1]
                     if m.sum() >= MIN_MUT else np.nan)
            d_lin = float(np.median(zz[m]) - np.median(zz[~m]))
            b = tmb.reindex(v.index).fillna(0).to_numpy()
            ry, rc, rb = (stats.rankdata(v.to_numpy()),
                          stats.rankdata(m.astype(float)), stats.rankdata(b))
            Z = np.column_stack([np.ones_like(rb), rb])
            ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
            ec = rc - Z @ np.linalg.lstsq(Z, rc, rcond=None)[0]
            pr = stats.pearsonr(ec, ey)
            rows.append({"label": label, "kind": kind, "allele": allele,
                         "drug": dn, "n_mut": int(m.sum()),
                         "n_wt": int((~m).sum()), "delta": d, "p": float(p),
                         "expect_sign": expect_sign,
                         "sign_ok": bool(np.sign(d) == expect_sign),
                         "p_lineage_adj": float(p_lin),
                         "sign_ok_lineage": bool(np.sign(d_lin) == expect_sign),
                         "p_tmb_adj": float(pr.pvalue),
                         "tmb_rho": float(pr.statistic), "note": ""})
        return rows

    all_rows = []
    print("\n=== POSITIVE CONTROLS (must pass for the negatives to mean "
          "anything) ===")
    for allele, drugs, sgn in CONTROLS:
        r = test(allele, drugs, sgn, allele, "control")
        all_rows += r
        for x in r:
            ok = ("PASS" if (x.get("sign_ok") and x["p"] < 0.05)
                  else ("wrong sign" if x.get("sign_ok") is False
                        else x["note"] or "n.s."))
            print(f"  {allele:18s} x {x['drug']:14s} n_mut={x['n_mut']:4d} "
                  f"delta={x['delta']:+.3f} p={x['p']:.2e}  {ok}"
                  if np.isfinite(x.get("p", np.nan)) else
                  f"  {allele:18s} x {x['drug']:14s} {x['note']}")

    # what does the class look like? build the EGFR-inhibitor list from GDSC
    tgt = (G.drop_duplicates("DRUG_NAME")
           .set_index("drug_l")[["PUTATIVE_TARGET", "PATHWAY_NAME"]])
    def class_drugs(patterns):
        """Match the drug's stated TARGET on word boundaries.

        Substring matching was tried first and is badly wrong here: "her"
        matches "other", so "Other, kinases" pulled 5-fluorouracil and
        5-azacytidine into the EGFR class and inflated it to 170 drugs. Targets
        are matched as whole tokens against PUTATIVE_TARGET only.
        """
        rx = re.compile(r"\b(" + "|".join(patterns) + r")\b", re.I)
        out = []
        for dn, r in tgt.iterrows():
            blob = str(r.PUTATIVE_TARGET)
            if blob.lower() in ("nan", "other", "unclassified"):
                continue
            if rx.search(blob):
                out.append(dn)
        return sorted(out)

    print("\n=== CANDIDATES ===")
    for cpd, allele, sgn, kw in CANDIDATES:
        lab = f"{allele} x {cpd}"
        exact = test(allele, [cpd], sgn, lab, "exact")
        all_rows += exact
        if exact and np.isfinite(exact[0].get("p", np.nan)):
            x = exact[0]
            verdict = ("REPLICATES" if x["p"] < 0.05 and x["sign_ok"]
                       else "does NOT replicate")
            print(f"  exact  {lab:34s} n_mut={x['n_mut']:3d} "
                  f"delta={x['delta']:+.3f} p={x['p']:.3f}  {verdict}")
        else:
            print(f"  exact  {lab:34s} compound absent from GDSC")
        cd = class_drugs(kw)
        if cd:
            cr = test(allele, cd, sgn, lab, "class")
            all_rows += cr
            ok = [x for x in cr if np.isfinite(x.get("p", np.nan))]
            if ok:
                sig = [x for x in ok if x["p"] < 0.05 and x["sign_ok"]]
                s_lin = [x for x in ok if x["p_lineage_adj"] < 0.05
                         and x["sign_ok_lineage"]]
                s_tmb = [x for x in ok if x["p_tmb_adj"] < 0.05]
                # Fisher assumes independence; drugs of one class tested on the
                # same lines are strongly correlated, so it is reported only as
                # a descriptive summary and the counts are what matter.
                print(f"  class  {len(ok)} {kw[0]}-class drugs tested: "
                      f"{len(sig)} raw, {len(s_lin)} after lineage adjustment, "
                      f"{len(s_tmb)} after TMB adjustment (expected by chance "
                      f"~{0.025*len(ok):.1f})")
                if sig:
                    print(f"         raw hits: "
                          f"{', '.join(x['drug'] for x in sig[:8])}")
    R = pd.DataFrame(all_rows)
    R.to_csv(TAB / "gdsc_allele_validation.csv", index=False)

    # Germline check. CCLE calls variants without matched normals, so common
    # polymorphisms survive into the somatic table. A germline SNP recurs at a
    # single genomic position across unrelated lines, carries no excess
    # mutational burden, and sits among other recurrent variants that are
    # SYNONYMOUS -- which no somatic hotspot does.
    print("\n=== GERMLINE CHECK on the surviving candidate alleles ===")
    for allele in sorted({a for _, a, _, _ in CANDIDATES}):
        gene = allele.split()[0]
        gm = dmg_all[dmg_all.Gene_symbol == gene]
        car = vsets.get(allele, set())
        if not car:
            continue
        top = gm.groupby("Protein_Change").RRID.nunique().sort_values(
            ascending=False).head(6)
        syn = [c for c in top.index
               if re.match(r"^p\.([A-Z])\d+\1$", str(c))]
        gc = gm[gm.Protein_Change.astype(str) ==
                allele.split(" ", 1)[1]].Genome_Change.nunique()
        b_car = tmb.reindex(list(car)).median()
        print(f"  {allele}: {len(car)} lines ({len(car)/tmb.shape[0]:.1%}), "
              f"{gc} distinct genomic position(s)")
        print(f"    mutational burden: carriers {b_car:.0f} vs all lines "
              f"{tmb.median():.0f}")
        print(f"    other most-recurrent {gene} variants: "
              f"{', '.join(str(c) for c in top.index[:5])}")
        if syn:
            print(f"    -> {len(syn)} of the top recurrent variants are "
                  f"SYNONYMOUS ({', '.join(syn[:3])}), the signature of "
                  f"germline\n       polymorphism surviving into an unmatched "
                  f"somatic call set. Treat as germline.")

    ctl = R[(R.kind == "control") & R.p.notna()]
    cand = R[(R.kind != "control") & R.p.notna()]
    n_ctl_pass = int(((ctl.p < 0.05) & ctl.sign_ok).sum())
    n_cand_pass = int(((cand.p < 0.05) & cand.sign_ok).sum())
    print(f"\nSUMMARY  controls: {n_ctl_pass}/{len(ctl)} replicate with the "
          f"correct sign")
    print(f"         candidates: {n_cand_pass}/{len(cand)} replicate with the "
          f"correct sign")
    if n_ctl_pass and not n_cand_pass:
        print("\n  The pipeline detects known biomarkers in GDSC and detects "
              "none of the\n  candidates. Combined with MUC16 and HMCN1 being "
              "the 2nd and 21st most\n  mutated genes of 18,739, the candidates "
              "should be treated as artefacts.")
    elif not n_ctl_pass:
        print("\n  Controls did not replicate — the validation is "
              "uninformative about the\n  candidates, and the mapping or metric "
              "should be checked before concluding.")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    plot = pd.concat([ctl.assign(grp="known biomarker"),
                      cand.assign(grp="novel candidate")])
    plot = plot[plot.p.notna()].copy()
    plot["lp"] = -np.log10(plot.p.clip(lower=1e-300))
    plot["signed"] = plot.lp * np.where(plot.sign_ok, 1, -1)
    yy = np.arange(len(plot))
    cols = [BLUE if g == "known biomarker" else ORANGE for g in plot.grp]
    ax[0].barh(yy, plot.signed, color=cols, height=0.7)
    ax[0].axvline(-np.log10(0.05), color="#444444", ls="--", lw=1)
    ax[0].set_yticks(yy, [f"{a.split()[0]} {a.split()[-1][:8]} × {d[:12]}"
                          for a, d in zip(plot.allele, plot.drug)], fontsize=6)
    ax[0].set_xlabel("−log10 p  (negative = wrong direction)")
    ax[0].set_title("A  Replication in GDSC", loc="left", fontweight="bold",
                    fontsize=10)

    ax[1].bar([0, 1], [n_ctl_pass / max(len(ctl), 1),
                       n_cand_pass / max(len(cand), 1)],
              color=[BLUE, ORANGE], width=0.55)
    ax[1].set_xticks([0, 1], [f"known biomarkers\n(n={len(ctl)})",
                              f"novel candidates\n(n={len(cand)})"])
    ax[1].set_ylabel("fraction replicating with correct sign")
    ax[1].set_ylim(0, 1)
    ax[1].set_title("B  The pipeline works; the candidates do not",
                    loc="left", fontweight="bold", fontsize=10)
    fig.suptitle("Do the novel allele candidates replicate in an independent "
                 "laboratory and assay?", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "gdsc_validation", FIG,
                    source_data={"validation": R}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
