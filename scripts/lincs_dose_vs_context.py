#!/usr/bin/env python3
"""Does the dose effect on context-dependence hold for a transcriptional readout?

PRISM shows it clearly for viability: context-dependence climbs with dose to a
peak near 2.5 uM and collapses at 10 uM, where every line is dying and the
shared component swamps everything. Tahoe cannot arbitrate -- its three
concentrations (0.05, 0.5, 5 uM) produce almost the same response magnitude
(shared variance 0.047 / 0.048 / 0.054), so there is no dose gradient to read.

LINCS phase 1 is the missing arm: transcriptional like Tahoe, but with a real
dose range and 2,834 compounds on replicate plates. If the inverted-U is a
property of the biology it should appear here; if it is a property of a killing
assay saturating, it should not.

The estimator is the same one used for PRISM and Tahoe -- covariance of per-line
residuals between independent replicate plates at matched compound and dose --
so the three are directly comparable.

Outputs: results/tables/dose_vs_context_lincs.csv
         figure bundle results/figures/13_prism/lincs_dose_vs_context/
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
LIN = ROOT / "data" / "external" / "lincs"
FIG = ROOT / "results" / "figures" / "13_prism"
TAB = ROOT / "results" / "tables"
MIN_LINES = 6
MIN_DOSES = 3
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gctx", default=str(LIN / "p1_level4.gctx"))
    ap.add_argument("--inst-info", default=str(LIN / "p1_inst_info.txt.gz"))
    ap.add_argument("--gene-info", default=str(LIN / "p1_gene_info.txt.gz"))
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    gi = pd.read_csv(args.gene_info, sep="\t")
    landmark = gi[gi.pr_is_lm == 1].pr_gene_id.astype(str).tolist()
    inst = pd.read_csv(args.inst_info, sep="\t", low_memory=False)
    inst["plate"] = inst.rna_plate.astype(str)
    cp = inst[(inst.pert_type == "trt_cp") & (inst.pert_dose > 0)].copy()
    # a compound is usable only if it spans several doses in several lines
    span = cp.groupby("pert_id").agg(nd=("pert_dose", "nunique"),
                                     nl=("cell_id", "nunique"))
    keep = span[(span.nd >= MIN_DOSES) & (span.nl >= MIN_LINES)].index
    cp = cp[cp.pert_id.isin(keep)]
    print(f"{len(cp)} instances, {cp.pert_id.nunique()} compounds spanning "
          f">={MIN_DOSES} doses in >={MIN_LINES} lines", flush=True)

    from cmapPy.pandasGEXpress.parse_gctx import parse
    g = parse(args.gctx, rid=landmark, cid=cp.inst_id.tolist())
    M = g.data_df.T
    meta = cp.set_index("inst_id").loc[M.index]
    X = M.to_numpy(dtype=np.float32)
    print(f"loaded {X.shape[0]} x {X.shape[1]}", flush=True)

    K = pd.DataFrame({"i": np.arange(len(meta)), "cpd": meta.pert_id.to_numpy(),
                      "line": meta.cell_id.to_numpy(),
                      "dose": meta.pert_dose.to_numpy(),
                      "plate": meta.plate.to_numpy()})
    # Each line's GENERAL transcriptional response to being perturbed at all,
    # taken leave-one-compound-out and removed before anything below is called
    # context-dependence. Like a genuine interaction it is shared between
    # replicate plates, so left in it enters the covariance and is counted as
    # context-specific. In PRISM the viability analogue reproduces across
    # disjoint compound halves at r = 0.989.
    # keyed by (line, PLATE): a Level 4 z-score is computed within its own
    # plate, so every profile a line has on a plate shares that plate's
    # normalisation. A correction pooled over plates would carry a share of each
    # and subtract it from the wrong side of the cross-plate pair below.
    _pc, _pk = {}, {}
    for (ln, pl, cp), g in K.groupby(["line", "plate", "cpd"], observed=True):
        _pc.setdefault((ln, pl), []).append(X[g.i.to_numpy()].mean(0))
        _pk.setdefault((ln, pl), []).append(cp)
    ALPHA = {}
    for key, vecs in _pc.items():
        V = np.stack(vecs)                       # compounds x genes
        tot, m = V.sum(0), len(V)
        if m >= 3:
            ALPHA[key] = {c: (tot - V[j]) / (m - 1)
                          for j, c in enumerate(_pk[key])}

    recs = []
    for cpd, gc in K.groupby("cpd", observed=True):
        doses = sorted(gc.dose.unique())
        if len(doses) < MIN_DOSES:
            continue
        for di, dose in enumerate(doses):
            gd = gc[gc.dose == dose]
            if gd.line.nunique() < MIN_LINES:
                continue
            # per-line mean profile at this dose, then the leave-one-line-out
            # shared response, exactly as in the PRISM and Tahoe versions
            lm, keys = [], []
            for ln, gl in gd.groupby("line", observed=True):
                lm.append(X[gl.i.to_numpy()].mean(0)); keys.append(ln)
            Lm = np.stack(lm)
            n = len(Lm)
            # subtract each line's general response, averaged over its plates
            # for the per-line mean profile (zero where not estimable)
            def _corr(ln, pl=None):
                if pl is not None:
                    return ALPHA.get((ln, pl), {}).get(cpd,
                                                       np.zeros(X.shape[1]))
                vs = [ALPHA[k][cpd] for k in ALPHA
                      if k[0] == ln and cpd in ALPHA[k]]
                return np.mean(vs, axis=0) if vs else np.zeros(X.shape[1])
            A = np.stack([_corr(ln) for ln in keys])
            A = A - A.mean(0, keepdims=True)
            Lm = Lm - A
            loo = (Lm.sum(0) - Lm) / (n - 1)
            add = float(np.mean(loo ** 2))
            cos = float(np.median([
                Lm[k] @ loo[k] / max(np.linalg.norm(Lm[k])
                                     * np.linalg.norm(loo[k]), 1e-9)
                for k in range(n)]))
            # interaction: residual agreement between replicate plates
            cov, npair = 0.0, 0
            for k, ln in enumerate(keys):
                gl = gd[gd.line == ln]
                pl = gl.plate.to_numpy(); ii = gl.i.to_numpy()
                for a in range(len(ii)):
                    for b in range(a + 1, len(ii)):
                        if pl[a] != pl[b]:
                            cov += float(np.mean(
                                (X[ii[a]] - loo[k] - _corr(ln, pl[a]))
                                * (X[ii[b]] - loo[k] - _corr(ln, pl[b]))))
                            npair += 1
            if npair >= 3:
                inter = max(cov / npair, 0.0)
                den = add + inter
                cdi = inter / den if den > 0 else np.nan
            else:
                inter, cdi = np.nan, np.nan
            recs.append({"compound": cpd, "dose": float(dose),
                         "dose_rank": di + 1, "n_dose": len(doses),
                         "n_lines": n, "n_pairs": npair, "additive": add,
                         "interaction": inter, "cdi": cdi, "cos_to_shared": cos})
    D = pd.DataFrame(recs)
    D["dose_pos"] = np.ceil(D.dose_rank / D.n_dose * 4).clip(1, 4).astype(int)
    D.to_csv(TAB / "dose_vs_context_lincs.csv", index=False)

    def trend(col):
        rs = []
        for _, g in D.groupby("compound", observed=True):
            g = g.dropna(subset=[col])
            if g.dose_rank.nunique() < 3:
                continue
            r = stats.spearmanr(g.dose_rank, g[col]).statistic
            if np.isfinite(r):
                rs.append(r)
        rs = np.array(rs)
        if len(rs) < 5:
            return None
        return {"n": len(rs), "rho": float(np.median(rs)),
                "pos": float((rs > 0).mean()),
                "p": float(stats.wilcoxon(rs).pvalue)}

    print(f"\n{D.compound.nunique()} compounds, "
          f"{D.cdi.notna().sum()} (compound, dose) cells with an estimable CDI")
    for col, lab in (("cdi", "context-dependence index"),
                     ("cos_to_shared", "alignment to the shared response"),
                     ("additive", "shared response magnitude")):
        t = trend(col)
        if t:
            print(f"  {lab:36s} vs dose: rho={t['rho']:+.3f}, "
                  f"{t['pos']:.0%} positive, p={t['p']:.2e} (n={t['n']})")
    S = D.groupby("dose_pos").agg(cdi=("cdi", "median"),
                                  add=("additive", "median"),
                                  inter=("interaction", "median"),
                                  cos=("cos_to_shared", "median"),
                                  dose=("dose", "median"),
                                  n=("compound", "size"))
    print(S.round(4).to_string())

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.4), constrained_layout=True)
    q = D.groupby("dose_pos").cdi.quantile([0.25, 0.5, 0.75]).unstack()
    ax[0].fill_between(q.index, q[0.25], q[0.75], color=VIOLET, alpha=0.2)
    ax[0].plot(q.index, q[0.5], "o-", color=VIOLET, lw=2, ms=7)
    ax[0].set_xticks(sorted(D.dose_pos.unique()))
    ax[0].set_xlabel("dose quartile (1 = lowest)")
    ax[0].set_ylabel("context-dependence index")
    t = trend("cdi")
    ax[0].set_title(f"A  LINCS phase 1 transcription\n"
                    f"rho={t['rho']:+.2f}, p={t['p']:.1e}" if t else "A",
                    loc="left", fontweight="bold", fontsize=10)

    ax[1].plot(S.index, S["add"], "o-", color=BLUE, lw=2, label="shared")
    ax[1].plot(S.index, S["inter"], "o-", color=ORANGE, lw=2,
               label="interaction")
    ax[1].set_yscale("log"); ax[1].set_xticks(sorted(D.dose_pos.unique()))
    ax[1].set_xlabel("dose quartile"); ax[1].set_ylabel("variance component")
    ax[1].legend(frameon=False, fontsize=8)
    ax[1].set_title("B  Which component grows faster", loc="left",
                    fontweight="bold", fontsize=10)

    qc = D.groupby("dose_pos").cos_to_shared.quantile([0.25, 0.5, 0.75]).unstack()
    ax[2].fill_between(qc.index, qc[0.25], qc[0.75], color=AQUA, alpha=0.2)
    ax[2].plot(qc.index, qc[0.5], "o-", color=AQUA, lw=2, ms=7)
    ax[2].set_xticks(sorted(D.dose_pos.unique()))
    ax[2].set_xlabel("dose quartile")
    ax[2].set_ylabel("alignment to the shared response")
    tc = trend("cos_to_shared")
    ax[2].set_title(f"C  Do lines converge or diverge?\n"
                    f"rho={tc['rho']:+.2f}, p={tc['p']:.1e}" if tc else "C",
                    loc="left", fontweight="bold", fontsize=10)
    fig.suptitle("Does dose change context-dependence for transcription too?",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "lincs_dose_vs_context", FIG,
                    source_data={"by_compound_dose": D,
                                 "summary": S.reset_index()}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
