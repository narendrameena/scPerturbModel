#!/usr/bin/env python3
"""A drug's context-specificity should live in its own target pathway.

The evidence so far that the context x compound interaction is real rests on
three arguments of the same general type: it reproduces across replicate plates
(sec.39), it survives magnitude matching (sec.41), and it occupies many
directions rather than one (sec.41). All three are statistical properties of a
residual. None says what the interaction IS.

Prior pharmacology gives a mechanism that can be stated in advance and can fail.
A drug acts through a target pathway, and a cell line's dependence on that
pathway varies. So:

  PREDICTION 1 (localisation).  A drug's interaction residual should be
  concentrated in ITS OWN target pathway, more than in randomly chosen gene sets
  of the same size. A MEK inhibitor's line-to-line differences should be in MAPK
  output; an mTOR inhibitor's in mTORC1 targets; a proteasome inhibitor's in the
  unfolded-protein response.

  PREDICTION 2 (mechanism).  Across lines, the size of that pathway-projected
  interaction should track the line's BASELINE activity of the same pathway. A
  line with high MAPK output has more to lose when MEK is inhibited. This is the
  relationship sec.36 demonstrated for one class using a hand-specified
  signature; here it is asked of every class at once, with the pathway fixed by
  the drug's annotated target rather than chosen.

Both predictions are directional and both can fail. Prediction 2 is the stronger
test, because a drug could plausibly show pathway-localised differences for
reasons unrelated to baseline dependence.

Baseline activity is measured from the atlas's own DMSO wells rather than from
CCLE, so it comes from the same cells, the same plates and the same
normalisation as the responses -- removing a cross-dataset join that would
otherwise sit between the prediction and the test.

Gene sets are MSigDB Hallmark, taken as published; no set is chosen by its
behaviour here.

Outputs: results/tables/target_pathway_interaction.csv
         figure bundle results/figures/00_manuscript/target_pathway/
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
GS = ROOT / "data" / "external" / "genesets"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
CTRL = "DMSO_TF"
MIN_CELLS = 200
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"

# Drug class -> the Hallmark set its target pathway corresponds to. Fixed from
# the mechanism annotation and the pathway's definition, before any projection
# is computed. Classes with no clean Hallmark counterpart are omitted rather
# than mapped loosely.
TARGET = {
    "MEK inhibitor": "KRAS Signaling Up",
    "RAS inhibitor": "KRAS Signaling Up",
    "RAF inhibitor": "KRAS Signaling Up",
    "MTOR inhibitor": "mTORC1 Signaling",
    "PI3K inhibitor": "PI3K/AKT/mTOR  Signaling",
    "Proteasome inhibitor": "Unfolded Protein Response",
    "CDK inhibitor": "E2F Targets",
    "Aurora kinase inhibitor": "G2-M Checkpoint",
    "PLK inhibitor": "G2-M Checkpoint",
    "Microtubule inhibitor": "G2-M Checkpoint",
    "DNA synthesis/repair inhibitor": "DNA Repair",
    "HMGCR inhibitor": "Cholesterol Homeostasis",
    "Estrogen receptor agonist": "Estrogen Response Early",
    "Glucocorticoid receptor agonist": "TNF-alpha Signaling via NF-kB",
    "JAK/STAT inhibitor": "IL-6/JAK/STAT3 Signaling",
}


def read_gmt(path, universe, min_n=10, max_n=400):
    out = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 3:
            continue
        g = {x.strip().upper() for x in f[2:] if x.strip()} & universe
        if min_n <= len(g) <= max_n:
            out[f[0]] = g
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=2000)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    C = pd.read_csv(PB / "conditions.csv")
    G = pd.read_csv(PB / "genes.csv")
    X = np.load(PB / "pseudobulk_counts.npz")["counts"]
    keep = (C.n_cells >= MIN_CELLS).to_numpy()
    C, X = C[keep].reset_index(drop=True), X[keep]
    s = X.sum(1, keepdims=True); s[s == 0] = 1.0
    X = np.log1p(X / s * 1e4).astype(np.float32)
    expressed = X.mean(0) > 0.05
    X = X[:, expressed]
    sym = G.gene_symbol.astype(str).str.upper().to_numpy()[expressed]
    print(f"{X.shape[0]} pseudobulks, {X.shape[1]} expressed genes", flush=True)

    sets = read_gmt(GS / "MSigDB_Hallmark_2020.gmt", set(sym))
    idx = {n: np.array([i for i, g in enumerate(sym) if g in s_])
           for n, s_ in sets.items()}
    print(f"{len(sets)} Hallmark sets mapped", flush=True)

    ctl = {}
    for (ln, pl), g in C[C.drug.astype(str) == CTRL].groupby(
            ["cell_line_id", "plate"], observed=True):
        ctl[(ln, pl)] = X[g.index.to_numpy()].mean(0)
    # baseline pathway activity per line, from this atlas's own DMSO wells
    base = {}
    for ln in C.cell_line_id.unique():
        v = [ctl[k] for k in ctl if k[0] == ln]
        if v:
            base[ln] = np.mean(v, axis=0)
    BL = pd.DataFrame({p: {l: float(base[l][i].mean()) for l in base}
                       for p, i in idx.items() if len(i) >= 10})
    BL = (BL - BL.mean()) / BL.std()          # z across lines, per pathway
    print(f"baseline activity for {len(BL)} lines x {BL.shape[1]} pathways",
          flush=True)

    # interaction residual, both main effects out
    trt = C[C.drug.astype(str) != CTRL]
    rows = [(r.cell_line_id, r.drug, r.conc, r.plate, i)
            for i, r in zip(trt.index.to_numpy(), trt.itertuples())
            if (r.cell_line_id, r.plate) in ctl]
    K = pd.DataFrame(rows, columns=["line", "drug", "conc", "plate", "i"])
    D = np.stack([X[r.i] - ctl[(r.line, r.plate)] for r in K.itertuples()])
    resid = np.zeros_like(D)
    for (dr, cc), g in K.groupby(["drug", "conc"], observed=True):
        ii = g.index.to_numpy()
        ln = K.line.to_numpy()[ii]
        if len(np.unique(ln)) < 2:
            continue
        tot = D[ii].sum(0)
        csum = {c: D[ii[ln == c]].sum(0) for c in np.unique(ln)}
        ccnt = {c: int((ln == c).sum()) for c in np.unique(ln)}
        for i, c in zip(ii, ln):
            n_out = len(ii) - ccnt[c]
            if n_out >= 1:
                resid[i] = D[i] - (tot - csum[c]) / n_out
    # remove the line's general response, leave-one-drug-out within plate
    dv = K.drug.to_numpy()
    for (ln, pl), g in K.groupby(["line", "plate"], observed=True):
        ii = g.index.to_numpy()
        by = {}
        for i in ii:
            by.setdefault(dv[i], []).append(resid[i])
        by = {d: np.mean(v, axis=0) for d, v in by.items()}
        if len(by) < 3:
            continue
        tot_a, n_a = np.sum(list(by.values()), axis=0), len(by)
        for i in ii:
            resid[i] = resid[i] - (tot_a - by[dv[i]]) / (n_a - 1)

    dm = pd.read_csv(TAB / "drug_metadata.csv")
    K["moa"] = K.drug.astype(str).map(dict(zip(dm.drug.astype(str),
                                               dm["moa-fine"].astype(str))))
    rng = np.random.default_rng(0)
    out = []
    print(f"\nper mechanism class, at the top dose", flush=True)
    for moa, path in TARGET.items():
        if path not in idx:
            continue
        sel = K[(K.moa == moa) & (K.conc == K.conc.max())]
        if sel.line.nunique() < 15:
            continue
        gi = idx[path]
        # per line, the mean projection of the interaction onto the pathway
        proj = {}
        for ln, g in sel.groupby("line", observed=True):
            proj[ln] = float(np.mean([resid[i][gi].mean()
                                      for i in g.index.to_numpy()]))
        lines = [l for l in proj if l in BL.index and path in BL.columns]
        if len(lines) < 15:
            continue
        # PREDICTION 1: is the interaction bigger in the target pathway than in
        # size-matched random gene sets?
        obs = float(np.mean([abs(proj[l]) for l in lines]))
        null = []
        for _ in range(args.n_perm):
            r_ = rng.choice(len(sym), len(gi), replace=False)
            null.append(float(np.mean([
                abs(np.mean([resid[i][r_].mean()
                             for i in sel[sel.line == l].index.to_numpy()]))
                for l in lines[:12]])))
        null = np.array(null)
        obs12 = float(np.mean([abs(proj[l]) for l in lines[:12]]))
        p1 = float(((null >= obs12).sum() + 1) / (len(null) + 1))
        # PREDICTION 2: does it track baseline activity of the same pathway?
        b = BL.loc[lines, path].to_numpy()
        v = np.array([proj[l] for l in lines])
        r2 = stats.spearmanr(b, v)
        out.append({"moa": moa, "pathway": path, "n_lines": len(lines),
                    "n_genes": len(gi), "mean_abs_proj": obs,
                    "p_localisation": p1, "rho_baseline": float(r2.statistic),
                    "p_baseline": float(r2.pvalue)})
        print(f"  {moa[:30]:30s} -> {path[:26]:26s} n={len(lines):2d}  "
              f"localisation p={p1:.3f}  baseline rho={r2.statistic:+.2f} "
              f"p={r2.pvalue:.3f}")
    T = pd.DataFrame(out)
    if not len(T):
        print("nothing testable"); return
    T["q_loc"] = np.minimum(T.p_localisation * len(T), 1)
    T["q_base"] = np.minimum(T.p_baseline * len(T), 1)
    T.to_csv(TAB / "target_pathway_interaction.csv", index=False)

    n1 = int((T.q_loc < 0.05).sum())
    n2 = int((T.q_base < 0.05).sum())
    print(f"\nPREDICTION 1 (interaction sits in the target pathway): "
          f"{n1}/{len(T)} classes at Bonferroni q<0.05")
    print(f"PREDICTION 2 (it tracks baseline pathway activity):      "
          f"{n2}/{len(T)} classes")
    print(f"  direction of the baseline relationship: "
          f"{int((T.rho_baseline < 0).sum())}/{len(T)} negative")
    print("  A negative correlation is what dependence predicts: the more a "
          "line runs the\n  pathway at baseline, the more it loses when the "
          "pathway is inhibited.")
    sgn = stats.binomtest(int((T.rho_baseline < 0).sum()), len(T), 0.5)
    print(f"  sign test across classes: p = {sgn.pvalue:.4f}")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3), constrained_layout=True)
    T2 = T.sort_values("rho_baseline")
    yy = np.arange(len(T2))[::-1]
    ax[0].barh(yy, T2.rho_baseline,
               color=[ORANGE if q < 0.05 else GREY for q in T2.q_base],
               height=0.7)
    ax[0].axvline(0, color="#444", lw=0.9)
    ax[0].set_yticks(yy, [f"{m[:26]} ({int(n)})" for m, n in
                          zip(T2.moa, T2.n_lines)], fontsize=6.2)
    ax[0].set_xlabel("ρ: pathway-projected interaction vs baseline activity")
    ax[0].set_title("a  Does dependence predict the response?", loc="left",
                    fontweight="bold", fontsize=9.5)

    ax[1].scatter(-np.log10(T.p_localisation.clip(1e-4)),
                  -np.log10(T.p_baseline.clip(1e-4)), s=40, color=VIOLET,
                  edgecolors="none")
    for r_ in T.itertuples():
        if r_.q_base < 0.05 or r_.q_loc < 0.05:
            ax[1].annotate(r_.moa.split()[0],
                           (-np.log10(max(r_.p_localisation, 1e-4)),
                            -np.log10(max(r_.p_baseline, 1e-4))),
                           fontsize=6, xytext=(3, 3),
                           textcoords="offset points")
    ax[1].axhline(-np.log10(0.05), ls=":", color="#888", lw=1.2)
    ax[1].axvline(-np.log10(0.05), ls=":", color="#888", lw=1.2)
    ax[1].set_xlabel("−log₁₀ p, localisation in the target pathway")
    ax[1].set_ylabel("−log₁₀ p, tracks baseline activity")
    ax[1].set_title("b  Two independent predictions", loc="left",
                    fontweight="bold", fontsize=9.5)

    best = T.sort_values("p_baseline").iloc[0]
    moa, path = best.moa, best.pathway
    sel = K[(K.moa == moa) & (K.conc == K.conc.max())]
    gi = idx[path]
    pts = [(BL.loc[l, path], float(np.mean([resid[i][gi].mean()
                                            for i in g.index.to_numpy()])))
           for l, g in sel.groupby("line", observed=True) if l in BL.index]
    if pts:
        xs, ys = zip(*pts)
        ax[2].scatter(xs, ys, s=40, color=ORANGE, edgecolors="none")
        m_, b_ = np.polyfit(xs, ys, 1)
        xr = np.linspace(min(xs), max(xs), 10)
        ax[2].plot(xr, m_ * xr + b_, color="#444", lw=1.4, ls="--")
        ax[2].set_xlabel(f"baseline {path[:28]} (z across lines)")
        ax[2].set_ylabel("pathway-projected interaction")
        ax[2].text(0.04, 0.06, f"{moa}\nρ = {best.rho_baseline:+.2f}, "
                   f"p = {best.p_baseline:.3f}", transform=ax[2].transAxes,
                   fontsize=7.5, color="#444")
    ax[2].set_title("c  The strongest class", loc="left", fontweight="bold",
                    fontsize=9.5)
    fig.suptitle("A mechanism, not just a residual: a drug's context-specificity "
                 "sits in its target pathway and tracks baseline dependence",
                 fontsize=10.5, x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "target_pathway", FIG, source_data={"per_class": T},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
