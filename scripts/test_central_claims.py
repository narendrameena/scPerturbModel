#!/usr/bin/env python3
"""The two central claims, as per-unit hypotheses with permutation nulls and FDR.

Both claims have so far been supported by pooled medians and a Wilcoxon over
compounds. That is the wrong inference twice over: the compounds are correlated
(effective n about 11 of 150), and a median tells you nothing about how many
individual units carry the effect. Each claim is restated here as a family of
hypotheses with its own null, and corrected across the family.

CLAIM 1 -- the variation is a drug x cell RELATION, not a cell property.

  If a cell line had a "response personality", its deviation from the compound
  average would be consistent across DIFFERENT compounds. So for each line we
  split the compounds into disjoint halves and correlate its mean residual
  between them. A line property predicts a positive correlation; a pure
  interaction predicts zero. The comparison quantity is the same line's
  interaction reproducibility across independent replicate plates at matched
  compound and dose, computed on the identical data.

  Null: permute the line labels within each compound, which destroys any
  line-level structure while preserving every marginal. FDR across lines.

CLAIM 2 -- it is governed by transcriptional STATE, not genotype.

  For each compound separately, is baseline expression predictive of the
  interaction residual, and is genotype? Cross-validated R^2 is compared against
  a null built by permuting the cell-line labels of y within that compound,
  which preserves the predictor block entirely and destroys only the pairing.
  Because ridge in the dual depends on y only through the right-hand side, the
  kernel is factorised once per fold and reused across permutations, which makes
  a per-compound permutation test affordable.

  FDR across compounds, separately for each predictor block. The headline is not
  a median but a COUNT: in how many compounds does each block beat its own null?

Outputs: results/tables/central_claim_line_vs_interaction.csv
         results/tables/central_claim_block_permutation.csv
         figure bundle results/figures/00_manuscript/central_claims/
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PR = ROOT / "data" / "external" / "prism"
FIG = ROOT / "results" / "figures" / "00_manuscript"
TAB = ROOT / "results" / "tables"
MIN_CPD = 40
N_PERM = 400
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def bh(p):
    p = np.asarray(p, float); n = len(p); q = np.empty(n); prev = 1.0
    for r, i in enumerate(np.argsort(p)[::-1]):
        prev = min(prev, p[i] * n / (n - r)); q[i] = prev
    return q


def cv_r2_perm(X, y, folds, alpha, n_perm, rng):
    """Observed CV R^2 and a permutation null, reusing one factorisation.

    Ridge in the dual is alpha = (K + lam I)^-1 y, so K depends only on X.
    Permuting y changes only the right-hand side, and the expensive solve can be
    reused -- which is what makes a per-compound permutation test tractable.
    """
    n = len(y)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    if ss_tot <= 0 or X.shape[1] == 0:
        return np.nan, np.nan, np.nan
    pre = []
    for f in np.unique(folds):
        tr, te = folds != f, folds == f
        A, B = X[tr], X[te]
        mu, sd = A.mean(0), A.std(0)
        sd = np.where(sd > 0, sd, 1.0)
        A = (A - mu) / sd; B = (B - mu) / sd
        Kk = A @ A.T + alpha * np.eye(A.shape[0])
        pre.append((tr, te, np.linalg.inv(Kk), B @ A.T))

    def score(yv):
        pred = np.zeros_like(yv)
        for tr, te, Kinv, BA in pre:
            ym = yv[tr].mean()
            pred[te] = ym + BA @ (Kinv @ (yv[tr] - ym))
        return 1 - float(np.sum((yv - pred) ** 2)) / ss_tot

    obs = score(y)
    null = np.array([score(y[rng.permutation(n)]) for _ in range(n_perm)])
    p = float((null.sum() * 0 + (null >= obs).sum() + 1) / (n_perm + 1))
    return obs, p, float(np.median(null))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-compounds", type=int, default=120)
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(ROOT / "scripts"))
    from genetic_architecture import residuals
    from expression_gap_closure import load_expression

    # Claim 1 is asked TWICE, of two residuals, because on the corrected one it
    # is circular: general sensitivity has been removed by construction, so
    # finding none proves only that the subtraction ran. The legacy residual is
    # where the question has content -- it is what every earlier analysis used.
    print("building PRISM residuals (legacy: compound mean only) ...",
          flush=True)
    resid, ti = residuals(legacy=True)

    # ---------------- CLAIM 1: line property vs interaction ----------------
    # assemble the (line, compound) residual matrix
    cpds = sorted(resid, key=lambda c: -len(resid[c]))[:600]
    M = pd.DataFrame({c: resid[c] for c in cpds})
    M = M.dropna(thresh=MIN_CPD)
    print(f"claim 1: {M.shape[0]} lines x {M.shape[1]} compounds", flush=True)

    rng = np.random.default_rng(0)
    cols = np.array(M.columns)
    perm = rng.permutation(len(cols))
    hA, hB = cols[perm[:len(cols) // 2]], cols[perm[len(cols) // 2:]]

    def line_consistency(mat, ha=hA, hb=hB):
        """Per line: does its mean residual on compound-half A agree with B?

        This is a single number per line; the test across lines asks whether
        those halves agree more than a line-label permutation allows. The
        halves are passed explicitly so the same split can be reused on a
        different residual matrix.
        """
        ca = [c for c in ha if c in mat.columns]
        cb = [c for c in hb if c in mat.columns]
        if not ca or not cb:
            return pd.Series(dtype=float), pd.Series(dtype=float)
        a = mat[ca].mean(axis=1)
        b = mat[cb].mean(axis=1)
        j = a.dropna().index.intersection(b.dropna().index)
        return a[j], b[j]

    a_obs, b_obs = line_consistency(M)
    r_obs = float(stats.spearmanr(a_obs, b_obs).statistic)
    # null: permute line labels independently within each compound
    null_r = []
    Mv = M.to_numpy()
    for _ in range(200):
        Z = Mv.copy()
        for k in range(Z.shape[1]):
            col = Z[:, k]
            ok = np.where(np.isfinite(col))[0]
            Z[ok, k] = col[rng.permutation(ok)]
        Zd = pd.DataFrame(Z, index=M.index, columns=M.columns)
        aa, bb = line_consistency(Zd)
        null_r.append(float(stats.spearmanr(aa, bb).statistic))
    null_r = np.array(null_r)
    p_line = float(((null_r >= r_obs).sum() + 1) / (len(null_r) + 1))
    print(f"\nCLAIM 1 -- is there a LINE-level component?")
    print(f"  agreement of a line's mean residual between disjoint compound "
          f"halves:\n    observed rho = {r_obs:+.3f}   null "
          f"{null_r.mean():+.3f} [{np.percentile(null_r,2.5):+.3f},"
          f"{np.percentile(null_r,97.5):+.3f}]   p = {p_line:.4f}")

    # per-line FDR: does THIS line have a reproducible personality?
    rows = []
    for ln in M.index:
        va = M.loc[ln, hA].dropna(); vb = M.loc[ln, hB].dropna()
        if len(va) < 15 or len(vb) < 15:
            continue
        # a line property means a non-zero mean shift consistent in both halves
        t1 = stats.ttest_1samp(va, 0.0); t2 = stats.ttest_1samp(vb, 0.0)
        same_sign = np.sign(va.mean()) == np.sign(vb.mean())
        rows.append({"line": ln, "mean_A": float(va.mean()),
                     "mean_B": float(vb.mean()), "same_sign": bool(same_sign),
                     "p": float(max(t1.pvalue, t2.pvalue))})
    L = pd.DataFrame(rows)
    if len(L):
        L["q"] = bh(L.p.to_numpy())
        L["line_property"] = (L.q < 0.05) & L.same_sign
        L.to_csv(TAB / "central_claim_line_vs_interaction.csv", index=False)
        print(f"  per-line test, FDR across {len(L)} lines: "
              f"{int(L.line_property.sum())} lines "
              f"({L.line_property.mean():.1%}) show a reproducible line-level "
              f"shift")
        print("  A line 'personality' would make this fraction large; an "
              "interaction-only\n  architecture predicts it near the FDR level.")

    # the same test on the CORRECTED residual, as a check that the subtraction
    # does what it claims and leaves no line-level term behind
    from perturbmodel.celldrug import prism_gamma, apportion
    gam, _, dec = prism_gamma()
    Mc = pd.DataFrame({c: gam[c] for c in cpds if c in gam}).dropna(
        thresh=MIN_CPD)
    ac, bc = line_consistency(Mc)
    r_corr = (float(stats.spearmanr(ac, bc).statistic) if len(ac) > 10
              else float("nan"))
    print(f"  same test on the CORRECTED residual: rho = {r_corr:+.3f} "
          f"(was {r_obs:+.3f})")

    # The claim "a relation, not a property" is settled by the SIZES of the two
    # terms, not by whether a property exists -- it plainly does. Both are
    # measured as covariances between independent estimates, so noise, which is
    # uncorrelated between them, contributes zero to either.
    v = apportion(dec)
    print(f"\n  variance apportionment, replicate-validated "
          f"({v['n_compounds']} compounds x {v['n_lines']} lines):")
    print(f"    drug effect                    {v['share_drug']:6.1%}")
    print(f"    cell property (general sens.)  {v['share_cell_property']:6.1%}"
          f"   var {v['var_cell_property']:.4f}")
    print(f"    cell-drug relation             {v['share_relation']:6.1%}"
          f"   var {v['var_cell_drug_relation']:.4f}")
    ratio = v['var_cell_drug_relation'] / max(v['var_cell_property'], 1e-9)
    print(f"  The relation is {ratio:.1f}x the property. Both are real; "
          f"'rather than a\n  cell property' overstates it, 'the larger of the "
          f"two' does not.")

    # ---------------- CLAIM 2: state vs genotype, per compound -------------
    EXPR = load_expression()
    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    dep = smp.cross_references.astype(str).str.extract(r"DepMap;\s*(ACH-\d+)")[0]
    cv2dep = dict(zip(smp.Cellosaurus_ID, dep))
    mut = pd.read_csv(PR / "mutation_long.tsv.gz", sep="\t", low_memory=False)
    mut["depmap"] = mut.RRID.map(cv2dep)
    mut = mut[mut.depmap.notna()]
    ns = mut[mut.Variant_Classification.astype(str).str.contains(
        "Missense|Nonsense|Frame_Shift|Splice|In_Frame", na=False)]
    gc = ns.groupby("Gene_symbol").depmap.nunique()
    genes = sorted(gc[gc >= 20].index)
    gsets = {g: set(d.depmap) for g, d in
             ns[ns.Gene_symbol.isin(genes)].groupby("Gene_symbol")}
    ci = pd.read_csv(PR / "secondary-screen-cell-line-info.csv")
    tis = dict(zip(ci.depmap_id, ci.primary_tissue.astype(str)))

    rows2 = []
    # Claim 2 is asked of the CORRECTED residual. On the legacy one a predictor
    # could win by predicting general sensitivity, which is a property of the
    # culture and says nothing about which drug it will respond to.
    sel = sorted(gam, key=lambda c: -len(gam[c]))[:args.n_compounds]
    for k, cpd in enumerate(sel):
        y0 = gam[cpd]
        keep = [d for d in y0.index if d in EXPR.index]
        if len(keep) < 80:
            continue
        yv = y0.loc[keep].to_numpy().astype(float)
        yv = yv - yv.mean()
        n = len(keep)
        folds = np.random.default_rng(k).permutation(np.arange(n) % 5)
        rr = np.random.default_rng(1000 + k)
        XE = EXPR.loc[keep].to_numpy(float)
        XE = XE[:, XE.std(0) > 0]
        XN = np.stack([np.isin(keep, list(gsets[g])).astype(float)
                       for g in genes], 1)
        XN = XN[:, XN.std(0) > 0]
        lin = np.array([tis.get(d, "?") for d in keep])
        ul = [u for u in np.unique(lin) if (lin == u).sum() >= 5]
        XL = np.stack([(lin == u).astype(float) for u in ul], 1) if ul \
            else np.zeros((n, 0))
        rec = {"compound": cpd, "n_lines": n}
        for lab, X in (("expression", XE), ("mutations", XN),
                       ("lineage", XL)):
            o, p, nm = cv_r2_perm(X, yv, folds, 1000.0, args.n_perm, rr)
            rec[f"r2_{lab}"] = o; rec[f"p_{lab}"] = p; rec[f"null_{lab}"] = nm
        rows2.append(rec)
        if k % 20 == 0:
            print(f"  {k+1}/{len(sel)} {cpd[:18]:18s} expr r2="
                  f"{rec['r2_expression']:+.3f} p={rec['p_expression']:.3f} | "
                  f"mut r2={rec['r2_mutations']:+.3f} "
                  f"p={rec['p_mutations']:.3f}", flush=True)
    B = pd.DataFrame(rows2)
    for lab in ("expression", "mutations", "lineage"):
        B[f"q_{lab}"] = bh(B[f"p_{lab}"].to_numpy())
    B.to_csv(TAB / "central_claim_block_permutation.csv", index=False)

    print(f"\nCLAIM 2 -- per-compound permutation test, FDR across "
          f"{len(B)} compounds:")
    for lab in ("expression", "lineage", "mutations"):
        sig = int((B[f"q_{lab}"] < 0.05).sum())
        print(f"  {lab:11s} {sig:3d}/{len(B)} compounds significant "
              f"({sig/len(B):5.1%})   median R2 {B[f'r2_{lab}'].median():+.4f}"
              f"   median null {B[f'null_{lab}'].median():+.4f}")
    ne, nm_ = int((B.q_expression < 0.05).sum()), int((B.q_mutations < 0.05).sum())
    tab = np.array([[int(((B.q_expression < 0.05) & (B.q_mutations < 0.05)).sum()),
                     int(((B.q_expression < 0.05) & (B.q_mutations >= 0.05)).sum())],
                    [int(((B.q_expression >= 0.05) & (B.q_mutations < 0.05)).sum()),
                     int(((B.q_expression >= 0.05) & (B.q_mutations >= 0.05)).sum())]])
    print(f"\n  expression vs mutations, per compound (rows: expr sig/not; "
          f"cols: mut sig/not):\n{tab}")
    if tab.sum() > 0:
        _, pm = stats.fisher_exact(tab) if tab.shape == (2, 2) else (0, np.nan)
        print(f"  McNemar-style contrast {ne} vs {nm_} compounds "
              f"(binomial p = "
              f"{stats.binomtest(tab[0,1], tab[0,1]+tab[1,0], 0.5).pvalue:.2e}"
              if (tab[0, 1] + tab[1, 0]) > 0 else "  (no discordant pairs)")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.1), constrained_layout=True)
    ax[0].hist(null_r, bins=30, color=GREY, alpha=0.85, label="line-label null")
    ax[0].axvline(r_obs, color=ORANGE, lw=2.4, label=f"observed {r_obs:+.2f}")
    ax[0].set_xlabel("agreement of a line's mean residual across compound halves")
    ax[0].set_ylabel("permutations"); ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_title("a  Is there a line-level component?", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(L):
        ax[1].bar([0, 1], [L.line_property.mean(), 0.05], width=0.55,
                  color=[VIOLET, GREY])
        ax[1].set_xticks([0, 1], [f"lines with a reproducible\nshift "
                                  f"(FDR<0.05)", "FDR level"], fontsize=7.5)
        ax[1].set_ylabel("fraction of lines")
        ax[1].set_title("b  Per-line test, FDR corrected", loc="left",
                        fontweight="bold", fontsize=9.5)

    labs = ["expression", "lineage", "mutations"]
    sig = [float((B[f"q_{l}"] < 0.05).mean()) for l in labs]
    ax[2].bar(range(3), sig, width=0.6, color=[VIOLET, AQUA, ORANGE])
    for i_, v in enumerate(sig):
        ax[2].text(i_, v + 0.015, f"{v:.0%}", ha="center", fontsize=9,
                   fontweight="bold")
    ax[2].axhline(0.05, ls="--", color="#888", lw=1.2, label="FDR level")
    ax[2].set_xticks(range(3), labs, fontsize=8)
    ax[2].set_ylabel("compounds beating their own permutation null")
    ax[2].legend(frameon=False, fontsize=7.5)
    ax[2].set_title(f"c  Per-compound, FDR across {len(B)}", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("The two central claims, tested per unit with permutation "
                 "nulls and FDR control", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    d = save_figure(fig, "central_claims", FIG,
                    source_data={"per_line": L if len(L) else pd.DataFrame(),
                                 "per_compound": B}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
