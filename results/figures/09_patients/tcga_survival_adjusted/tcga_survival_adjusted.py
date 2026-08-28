#!/usr/bin/env python3
"""Does the survival association of the response programs survive adjustment?

The unadjusted analysis found C1 predicts worse and C7 better overall survival
across cancer types. The first objection is confounding by disease burden and
sample composition, so this refits with covariates:

  age, sex, tumour stage (ordinal I-IV), histological grade (ordinal G1-G4),
  and tumour-purity proxies.

Purity note: the TCGA ABSOLUTE purity table is not publicly downloadable from
the Xena/GDC endpoints available here, so purity is proxied by transcriptome
based stromal and immune infiltration scores (ESTIMATE-style: mean z of
hallmark EMT genes and of hallmark immune genes). This adjusts for the dominant
axis of non-tumour content but is weaker than ABSOLUTE and is reported as such.

Model: Cox proportional hazards fitted separately WITHIN each cancer type
(so lineage is fully stratified), then the component's coefficient is combined
across cancer types by Stouffer, weighting by sqrt(events).

Outputs: results/tables/tcga_survival_adjusted.csv
         figure bundle results/figures/09_patients/tcga_survival_adjusted/
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
TCGA = ROOT / "data" / "external" / "tcga"
GS = ROOT / "data" / "external" / "genesets"
FIG = ROOT / "results" / "figures" / "09_patients"
TAB = ROOT / "results" / "tables"
MIN_EVENTS = 20
MIN_N = 60
BLUE, ORANGE = "#2a78d6", "#eb6834"

STAGE_MAP = {"I": 1, "II": 2, "III": 3, "IV": 4}
GRADE_MAP = {"G1": 1, "G2": 2, "G3": 3, "G4": 4}


def parse_stage(s):
    if not isinstance(s, str):
        return np.nan
    s = s.upper().replace("STAGE", "").strip()
    for k in ("IV", "III", "II", "I"):
        if s.startswith(k):
            return STAGE_MAP[k]
    return np.nan


def parse_grade(s):
    if not isinstance(s, str):
        return np.nan
    return GRADE_MAP.get(s.strip().upper()[:2], np.nan)


def hallmark_genes(name_contains, universe):
    out = set()
    for line in (GS / "MSigDB_Hallmark_2020.gmt").read_text().splitlines():
        p = line.rstrip().split("\t")
        if p and name_contains.lower() in p[0].lower():
            out |= {g.split(",")[0].strip().upper() for g in p[2:] if g.strip()}
    return sorted(out & universe)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="dev47")
    args = ap.parse_args()
    from lifelines import CoxPHFitter

    S = pd.read_csv(TAB / "tcga_component_scores.csv", index_col=0)
    comps = [c for c in S.columns if c.startswith("C") and c != "cancer_type"]
    surv = pd.read_csv(TCGA / "survival.tsv", sep="\t", index_col=0)
    surv = surv.loc[[i for i in S.index if i in surv.index]]
    S = S.loc[surv.index]

    print("loading expression for purity proxies ...", flush=True)
    expr = pd.read_csv(TCGA / "pancan_expr.tsv.gz", sep="\t", index_col=0)
    expr.index = expr.index.str.upper()
    expr = expr[~expr.index.duplicated()]
    expr = expr[[c for c in S.index if c in expr.columns]]
    uni = set(expr.index)
    strom = hallmark_genes("Epithelial Mesenchymal Transition", uni)
    imm = hallmark_genes("Inflammatory Response", uni)
    zz = expr.sub(expr.mean(axis=1), axis=0).div(
        expr.std(axis=1).replace(0, np.nan), axis=0)
    purity = pd.DataFrame({"stromal": zz.loc[strom].mean(),
                           "immune": zz.loc[imm].mean()})
    print(f"purity proxies from {len(strom)} stromal and {len(imm)} immune genes")
    del expr, zz

    df = S.join(purity)
    df["age"] = pd.to_numeric(surv["age_at_initial_pathologic_diagnosis"],
                              errors="coerce")
    df["male"] = (surv["gender"].astype(str).str.upper() == "MALE").astype(float)
    df["stage"] = surv["ajcc_pathologic_tumor_stage"].map(parse_stage)
    df["grade"] = surv["histological_grade"].map(parse_grade)
    df["time"] = pd.to_numeric(surv["OS.time"], errors="coerce")
    df["event"] = pd.to_numeric(surv["OS"], errors="coerce")
    df["ctype"] = S["cancer_type"]

    recs = []
    for comp in comps:
        for model, covs in (("unadjusted", []),
                            ("adjusted", ["age", "male", "stage", "grade",
                                          "stromal", "immune"])):
            zs, ws, used = [], [], 0
            for ct, sub in df.groupby("ctype"):
                cols = [comp, "time", "event"] + covs
                d = sub[cols].replace([np.inf, -np.inf], np.nan).dropna()
                d = d[d.time > 0]
                # drop covariates that are constant within this cancer type
                keep = [c for c in cols
                        if c in ("time", "event") or d[c].nunique() > 1]
                d = d[keep]
                ev = int(d.event.sum())
                if len(d) < MIN_N or ev < MIN_EVENTS or comp not in d.columns:
                    continue
                try:
                    cph = CoxPHFitter(penalizer=0.1)
                    cph.fit(d, duration_col="time", event_col="event")
                    z = float(cph.summary.loc[comp, "z"])
                except Exception:
                    continue
                zs.append(z); ws.append(np.sqrt(ev)); used += 1
            if used >= 5:
                w = np.array(ws)
                z = float(np.sum(w * np.array(zs)) / np.sqrt(np.sum(w ** 2)))
                recs.append({"component": comp, "model": model,
                             "n_cancer_types": used, "z": z,
                             "p": float(2 * stats.norm.sf(abs(z)))})
    res = pd.DataFrame(recs)
    if not len(res):
        print("no models fitted")
        return
    piv = res.pivot(index="component", columns="model", values="z")
    pp = res.pivot(index="component", columns="model", values="p")
    out = piv.join(pp, lsuffix="_z", rsuffix="_p")
    out.to_csv(TAB / "tcga_survival_adjusted.csv")
    print("\nCox z per component (positive = high score, worse survival):")
    print(out.round(4).to_string())

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(figsize=(8.5, 4.4), constrained_layout=True)
    x = np.arange(len(piv))
    ax.bar(x - 0.2, piv["unadjusted"], width=0.4, color="#bdbdbd",
           label="unadjusted")
    ax.bar(x + 0.2, piv["adjusted"], width=0.4, color=BLUE,
           label="adjusted (age, sex, stage, grade, stromal, immune)")
    for v in (-1.96, 1.96):
        ax.axhline(v, color="#888888", ls="--", lw=0.9)
    ax.set_xticks(x, piv.index)
    ax.set_ylabel("Cox z, combined across cancer types")
    ax.set_title("Survival association before and after adjustment",
                 loc="left", fontweight="bold", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    d = save_figure(fig, "tcga_survival_adjusted", FIG, source_data=res,
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
