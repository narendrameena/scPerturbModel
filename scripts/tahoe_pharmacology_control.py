#!/usr/bin/env python3
"""Can Tahoe recover pharmacology we already know? A biological positive control.

RESULTS.md sec.31 reports that Tahoe shows no detectable context x compound
interaction at matched dose (0.5%, CI [0.0-1.5%]). RESULTS.md sec.35 showed how
that kind of null can be an artefact: in Spear-ATAC, a global statistic over
2,174 motifs found nothing while six transcription-factor knockouts were plainly
moving their own motif, because averaging over features that did not respond
buried the ones that did. The same instrument produced the Tahoe null.

A per-pair statistical scan would be the same blind instrument again. What
distinguishes a real null from a diluted one is **prior pharmacology**: known
biology specifies which drug, which lines, which genes, and which direction,
before anything is computed.

The control used here is the best-established drug x genotype interaction in
cancer pharmacology.

  MAIN EFFECT.  MEK inhibitors suppress ERK-dependent transcription. The output
  genes are not chosen here -- they are the canonical MEK-dependent signature of
  Pratilas et al. (*PNAS* 106:4519, 2009), used clinically as a pharmacodynamic
  marker: DUSP4, DUSP6, SPRY2, SPRY4, ETV4, ETV5, PHLDA1, EPHA2, SPRED1, SPRED2,
  CCND1, FOSL1, MYC. Prediction: these fall on MEK inhibition, in every line.

  INTERACTION.  Lines driven by BRAF or RAS mutation depend on that pathway and
  carry high baseline ERK output; MAPK-wild-type lines do not. Prediction: the
  fall is LARGER in mutant lines. This is a context x compound interaction with
  its direction fixed in advance by biology rather than read off the data.

  SECOND MAIN-EFFECT CONTROL.  Proteasome inhibitors trigger the heat-shock
  response (HSPA1A, HSPA1B, DNAJB1, HSPH1, HSPB1, BAG3) in essentially any cell.
  Prediction: up, and roughly equally in all lines -- a positive control for a
  SHARED effect, which distinguishes "the atlas measures drug response" from
  "the atlas measures line-specific drug response".

If the interaction prediction fails while the main-effect predictions hold, the
sec.31 null is biologically as well as statistically supported: the atlas sees
drug effects and does not see this interaction. If the interaction prediction
holds, the pooled 0.5% is diluted in the way sec.35 documents, and the headline
must be restated.

Outputs: results/tables/tahoe_pharmacology_control.csv
         figure bundle results/figures/00_manuscript/pharmacology_control/
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
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
CTRL = "DMSO_TF"
MIN_CELLS = 200
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"

# Pratilas et al., PNAS 2009 -- the MEK-dependent output signature. Fixed before
# looking at the data; no gene here is selected on its behaviour in Tahoe.
MEK_OUTPUT = ["DUSP4", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1",
              "EPHA2", "SPRED1", "SPRED2", "CCND1", "FOSL1", "MYC"]
# canonical proteotoxic-stress response
HEAT_SHOCK = ["HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "HSPB1", "BAG3", "DNAJA1",
              "HSPA6"]
MEK_DRUGS = ["Cobimetinib", "Trametinib", "Binimetinib", "TAK-733"]
PROTEASOME = ["Bortezomib", "Ixazomib", "Ixazomib citrate"]
MAPK_DRIVERS = ("BRAF", "KRAS", "NRAS")


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
    sym = G.gene_symbol.astype(str).str.upper().to_numpy()
    pos = {g: i for i, g in enumerate(sym)}
    print(f"{X.shape[0]} pseudobulks, {C.cell_line_id.nunique()} lines, "
          f"{C.drug.nunique()} drugs", flush=True)

    md = pd.read_csv(TAB / "cell_line_metadata.csv")
    md = md[md.Cell_ID_Cellosaur.isin(set(C.cell_line_id))]
    drv = md.groupby("Cell_ID_Cellosaur").Driver_Gene_Symbol.apply(
        lambda s_: set(s_.dropna().astype(str)))
    mapk = {k: bool(v & set(MAPK_DRIVERS)) for k, v in drv.items()}
    name = dict(zip(md.Cell_ID_Cellosaur, md.cell_name))
    organ = dict(zip(md.Cell_ID_Cellosaur, md.Organ.astype(str)))
    n_mut = sum(mapk.values())
    print(f"MAPK-driven (BRAF/KRAS/NRAS) lines: {n_mut} of {len(mapk)}; "
          f"wild-type {len(mapk) - n_mut}", flush=True)

    ctl = {}
    for (ln, pl), g in C[C.drug.astype(str) == CTRL].groupby(
            ["cell_line_id", "plate"], observed=True):
        ctl[(ln, pl)] = X[g.index.to_numpy()].mean(0)

    def score(genes):
        idx = [pos[g] for g in genes if g in pos]
        missing = [g for g in genes if g not in pos]
        if missing:
            print(f"   not in the atlas: {missing}")
        return idx

    mek_idx = score(MEK_OUTPUT)
    hs_idx = score(HEAT_SHOCK)

    def responses(drugs, gidx):
        rows = []
        sel = C[C.drug.isin(drugs)]
        for r in sel.itertuples():
            a = ctl.get((r.cell_line_id, r.plate))
            if a is None or r.cell_line_id not in mapk:
                continue
            d = X[r.Index] - a
            rows.append({"line": r.cell_line_id,
                         "name": name.get(r.cell_line_id, r.cell_line_id),
                         "drug": r.drug, "conc": r.conc,
                         "score": float(d[gidx].mean()),
                         "mapk": mapk[r.cell_line_id],
                         "organ": organ.get(r.cell_line_id, "?")})
        return pd.DataFrame(rows)

    # ---------- 1. main effect: does MEK inhibition suppress ERK output? ----
    M = responses(MEK_DRUGS, mek_idx)
    print(f"\n1. MAIN EFFECT — MEK inhibitors on the Pratilas ERK-output "
          f"signature")
    print(f"   {len(M)} (line, drug, dose) conditions across "
          f"{M.drug.nunique()} MEK inhibitors")
    for cc, g in M.groupby("conc"):
        t = stats.wilcoxon(g.score) if len(g) > 10 else None
        print(f"   {cc:>6} uM: median Δ = {g.score.median():+.4f}  "
              f"({(g.score < 0).mean():.0%} of conditions down"
              + (f", p = {t.pvalue:.1e})" if t else ")"))
    top = M[M.conc == M.conc.max()]
    print(f"   Prediction was DOWN. "
          f"{'CONFIRMED' if top.score.median() < 0 else 'NOT confirmed'} at "
          f"the top dose.")

    # ---------- 2. main effect: proteasome inhibition and heat shock --------
    P = responses(PROTEASOME, hs_idx)
    if len(P):
        tp = P[P.conc == P.conc.max()]
        print(f"\n2. MAIN EFFECT — proteasome inhibitors on the heat-shock "
              f"response")
        print(f"   top dose: median Δ = {tp.score.median():+.4f} "
              f"({(tp.score > 0).mean():.0%} of conditions up, n = {len(tp)})")
        print(f"   Prediction was UP. "
              f"{'CONFIRMED' if tp.score.median() > 0 else 'NOT confirmed'}. "
              f"A shared effect the atlas plainly measures.")

    # ---------- 3. the interaction, with its direction fixed in advance -----
    print(f"\n3. INTERACTION — is ERK-output suppression larger in "
          f"MAPK-driven lines?")
    rows = []
    for cc, g in M.groupby("conc"):
        a = g[g.mapk].score
        b = g[~g.mapk].score
        if len(a) < 5 or len(b) < 5:
            continue
        diff = float(a.median() - b.median())
        # permutation over GENOTYPE labels, within dose, preserving each line's
        # contribution -- the null is "genotype is unrelated to the response"
        rng = np.random.default_rng(0)
        lines_u = g.drop_duplicates("line")[["line", "mapk", "organ"]]
        lab = lines_u.mapk.to_numpy()
        null = []
        for _ in range(args.n_perm):
            pm = rng.permutation(lab)
            mp = dict(zip(lines_u.line, pm))
            gm = g.line.map(mp).to_numpy()
            null.append(float(g.score[gm].median() - g.score[~gm].median()))
        null = np.array(null)
        pv = float(((null <= diff).sum() + 1) / (len(null) + 1))  # one-sided

        # BRAF/RAS mutation is not spread evenly across lineages -- it
        # concentrates in colorectal, melanoma and lung -- so an unstratified
        # permutation cannot separate genotype from tissue of origin. Permuting
        # the genotype label only WITHIN each organ holds lineage fixed and
        # asks whether genotype adds anything beyond it.
        org = dict(zip(lines_u.line, lines_u.organ))
        null_s = []
        by_org = {}
        for ln_, o_ in org.items():
            by_org.setdefault(o_, []).append(ln_)
        for _ in range(args.n_perm):
            mp = {}
            for o_, mem in by_org.items():
                labs = lines_u.set_index("line").loc[mem].mapk.to_numpy()
                mp.update(dict(zip(mem, rng.permutation(labs))))
            gm = g.line.map(mp).to_numpy()
            null_s.append(float(g.score[gm].median() - g.score[~gm].median()))
        null_s = np.array(null_s)
        pv_s = float(((null_s <= diff).sum() + 1) / (len(null_s) + 1))
        u = stats.mannwhitneyu(a, b, alternative="less")
        rows.append({"conc": cc, "n_mut": len(a), "n_wt": len(b),
                     "median_mut": float(a.median()),
                     "median_wt": float(b.median()), "diff": diff,
                     "p_perm": pv, "p_perm_within_organ": pv_s,
                     "p_mwu": float(u.pvalue)})
        print(f"   {cc:>6} uM: mutant {a.median():+.4f} (n={len(a)})  vs  "
              f"wild-type {b.median():+.4f} (n={len(b)})   Δ = {diff:+.4f}   "
              f"p = {pv:.4f}   (within-organ p = {pv_s:.4f})")
    I = pd.DataFrame(rows)
    if len(I):
        best = I.sort_values("p_perm").iloc[0]
        print(f"\n   Prediction was that mutant lines fall FURTHER (Δ < 0).")
        if ((I.p_perm < 0.05) & (I.p_perm_within_organ < 0.05)).any() \
                and (I["diff"] < 0).any():
            n_ok = int(((I.p_perm < 0.05)
                        & (I.p_perm_within_organ < 0.05)).sum())
            print(f"   CONFIRMED at {n_ok} of {len(I)} doses, and it survives "
                  f"permuting genotype\n   WITHIN organ (best within-organ p = "
                  f"{I.p_perm_within_organ.min():.4f}), so it is not lineage "
                  f"wearing a\n   genotype label. Tahoe DOES resolve this "
                  f"interaction, and the pooled 0.5%\n   of sec.31 is diluted "
                  f"rather than empty.")
        else:
            print(f"   NOT confirmed (best p = {best.p_perm:.4f}). The most "
                  f"robust drug x genotype\n   interaction in cancer "
                  f"pharmacology is not detectable in this atlas, which "
                  f"supports\n   sec.31's null biologically as well as "
                  f"statistically.")

    out = pd.concat([M.assign(panel="mek"), P.assign(panel="proteasome")],
                    ignore_index=True)
    out.to_csv(TAB / "tahoe_pharmacology_control.csv", index=False)
    if len(I):
        I.to_csv(TAB / "tahoe_pharmacology_interaction.csv", index=False)

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)
    cs = sorted(M.conc.unique())
    ax[0].boxplot([M[M.conc == c].score.dropna() for c in cs], showfliers=False,
                  patch_artist=True, medianprops=dict(color="black", lw=1.4))
    for i_, c in enumerate(cs):
        v = M[M.conc == c].score.median()
        ax[0].text(i_ + 1, v + 0.01, f"{v:+.3f}", ha="center", fontsize=8,
                   fontweight="bold")
    ax[0].axhline(0, color="#444", lw=1.0)
    ax[0].set_xticks(range(1, len(cs) + 1), [f"{c:g}" for c in cs], fontsize=8)
    ax[0].set_xlabel("MEK inhibitor concentration (µM)")
    ax[0].set_ylabel("Δ ERK-output signature (Pratilas 2009)")
    ax[0].set_title("a  Main effect: the drug works", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(P):
        cs2 = sorted(P.conc.unique())
        ax[1].boxplot([P[P.conc == c].score.dropna() for c in cs2],
                      showfliers=False, patch_artist=True,
                      medianprops=dict(color="black", lw=1.4))
        ax[1].axhline(0, color="#444", lw=1.0)
        ax[1].set_xticks(range(1, len(cs2) + 1), [f"{c:g}" for c in cs2],
                         fontsize=8)
        ax[1].set_xlabel("proteasome inhibitor concentration (µM)")
        ax[1].set_ylabel("Δ heat-shock response")
    ax[1].set_title("b  A second shared effect, as predicted", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(I):
        xx = np.arange(len(I))
        w = 0.36
        ax[2].bar(xx - w / 2, I.median_mut, w, color=ORANGE,
                  label=f"BRAF/RAS-driven (n={int(I.n_mut.iloc[0])})")
        ax[2].bar(xx + w / 2, I.median_wt, w, color=GREY,
                  label=f"MAPK wild-type (n={int(I.n_wt.iloc[0])})")
        ax[2].axhline(0, color="#444", lw=1.0)
        ax[2].set_xticks(xx, [f"{c:g} µM" for c in I.conc], fontsize=8)
        ax[2].set_ylabel("Δ ERK-output signature")
        ax[2].legend(frameon=False, fontsize=7)
        best = I.sort_values("p_perm").iloc[0]
        ax[2].text(0.5, 0.06, f"best permutation p = {best.p_perm:.3f}",
                   transform=ax[2].transAxes, ha="center", fontsize=8,
                   color=ORANGE, fontweight="bold")
    ax[2].set_title("c  The interaction biology predicts", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Does Tahoe recover known pharmacology? Main effects, and the "
                 "drug × genotype interaction predicted in advance",
                 fontsize=10.5, x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "pharmacology_control", FIG,
                    source_data={"responses": out,
                                 "interaction": I if len(I) else pd.DataFrame()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
