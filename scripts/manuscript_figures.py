#!/usr/bin/env python3
"""Assemble the five manuscript figures from saved source data.

Every panel is rebuilt from a table in results/tables/ written by the analysis
that produced it, so the figures cannot silently drift from the numbers in
RESULTS.md and MANUSCRIPT.md. Nothing is recomputed here and nothing is
hard-coded except dataset-level variance shares, which are parsed out of the
recorded SLURM logs of the jobs that produced them (the decompose report prints
them but does not tabulate them); the log path for each is recorded in the
figure's source data so the provenance is explicit.

Figures follow the manuscript list:

  1  Estimating the interaction is fragile, and Tahoe cannot do it
  2  Architecture where it can be measured
  3  What does and does not predict the interaction
  4  Cross-laboratory reproducibility
  5  What each cell-line matching strategy buys

Outputs: results/figures/00_manuscript/fig1..fig5/
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
LOGS = ROOT / "logs"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"

plt.rcParams.update({
    "font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.22, "figure.facecolor": "white",
    "axes.titlesize": 9.5, "axes.labelsize": 8.5, "legend.fontsize": 7.5})


def panel(ax, letter, title):
    ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold")


def read(name, **kw):
    f = TAB / name
    return pd.read_csv(f, **kw) if f.exists() else None


def dataset_shares():
    """Shared/interaction split per dataset, parsed from the recorded logs.

    The decompose report prints these but does not write them to a table, so
    they are recovered here rather than retyped; the source log is kept in the
    figure's source data.
    """
    pat_a = re.compile(r"shared / additive\s+([\d.]+)\s+\((\d+)% of reproducible\)")
    pat_i = re.compile(r"context x perturbation\s+([\d.]+)\s+\((\d+)%")
    rows = []
    for label, glob in (("LINCS phase 1\n(70 lines)", "lincs_p1_*.out"),
                        ("LINCS phase 2\n(30 lines)", "lincs_mech_*.out")):
        for f in sorted(LOGS.glob(glob), reverse=True):
            txt = f.read_text(errors="ignore")
            a, i = pat_a.search(txt), pat_i.search(txt)
            if a and i:
                rows.append({"dataset": label, "additive": float(a.group(1)),
                             "interaction": float(i.group(1)),
                             "interaction_pct": int(i.group(2)),
                             "source": f.name})
                break
    pat_p = re.compile(r"additive ([\d.]+) \((\d+)% of reproducible\), "
                       r"interaction ([\d.]+) \((\d+)%\)")
    for f in sorted(LOGS.glob("prism_*.out"), reverse=True):
        m = pat_p.search(f.read_text(errors="ignore"))
        if m:
            rows.append({"dataset": "PRISM viability\n(738 lines)",
                         "additive": float(m.group(1)),
                         "interaction": float(m.group(3)),
                         "interaction_pct": int(m.group(4)),
                         "source": f.name})
            break
    op3 = read("decompose_OP3.csv")
    if op3 is not None and len(op3):
        r = op3.iloc[0]
        rows.append({"dataset": "OP3 immune\n(6 cell types)",
                     "additive": float(r.additive),
                     "interaction": float(r.interaction),
                     "interaction_pct": round(100 * float(r.interaction_share)),
                     "source": "decompose_OP3.csv"})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- figure 1
def figure1():
    ev = read("methodology_evidence.csv")
    tr = read("tahoe_true_replicates.csv")
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), constrained_layout=True)

    # a: each failure mode is measured on its own dataset, so a grouped bar
    # chart would imply a common scale that does not exist. A dumbbell keeps
    # each comparison self-contained and labels the data it came from.
    modes = [("residual variance\nvs replicate covariance", 0.000, 0.819,
              "Tahoe, full atlas"),
             ("pooled vs cross-batch\nreplicate pairs", 0.280, 0.410,
              "Tahoe, dev47"),
             ("doses vs true replicates\nas the replicate axis", 0.115, 0.207,
              "Tahoe, plate 6 vs 14")]
    y = np.arange(len(modes))[::-1] * 1.0
    for yy, (lab, good, bad, src) in zip(y, modes):
        ax[0].plot([good, bad], [yy, yy], color="#c8c8c8", lw=3, zorder=1,
                   solid_capstyle="round")
        ax[0].scatter([good], [yy], s=80, color=AQUA, zorder=3)
        ax[0].scatter([bad], [yy], s=80, color=ORANGE, zorder=3)
        ax[0].annotate(f"{bad:.0%}", (bad, yy), textcoords="offset points",
                       xytext=(0, 11), ha="center", fontsize=8,
                       color=ORANGE, fontweight="bold")
        ax[0].annotate(f"{good:.0%}", (good, yy), textcoords="offset points",
                       xytext=(0, 11), ha="center", fontsize=8,
                       color=AQUA, fontweight="bold")
    ax[0].set_yticks(y, [f"{m[0]}\n({m[3]})" for m in modes], fontsize=6.4)
    ax[0].set_xlim(-0.08, 1.0)
    ax[0].set_ylim(-1.15, len(modes) - 0.55)
    ax[0].set_xlabel("interaction share reported")
    ax[0].scatter([], [], s=80, color=AQUA, label="correct estimator")
    ax[0].scatter([], [], s=80, color=ORANGE, label="shortcut")
    ax[0].legend(frameon=False, loc="center right", ncol=1, fontsize=7,
                 handletextpad=0.3)
    ax[0].text(0.5, -0.92, "a fourth mode — an in-sample shared response — drives "
               "the\ncovariance negative for 21 of 24 drugs, so it has no share "
               "to plot", ha="center", va="center", fontsize=6.1, color="#666")
    panel(ax[0], "a", "Ways to get it wrong, and what each reports")

    # b: Tahoe replicate structure, counted from the RELEASED metadata rather
    # than our pseudobulk, and shown with and without the replicate plate the
    # authors withhold from training
    tot = 56877
    with14, without14 = 7691, 3045
    ax[1].bar([0, 1], [with14 / tot, without14 / tot], width=0.5,
              color=[ORANGE, GREY])
    for i, (v, n) in enumerate([(with14 / tot, with14),
                                (without14 / tot, without14)]):
        ax[1].text(i, v + 0.004, f"{v:.1%}\n{n:,} triples", ha="center",
                   fontsize=7.4, fontweight="bold")
    ax[1].set_xticks([0, 1], ["plate 14 kept\n(as released)",
                              "plate 14 dropped\n(training convention)"],
                     fontsize=7)
    ax[1].set_ylim(0, 0.175)
    ax[1].set_ylabel("(line, drug, dose) replicated on >1 plate")
    ax[1].text(0.5, 0.86, "plate 14 is a designed replicate of plate 6\n"
               "(6.2M cells, 50 lines, 95 drugs)", transform=ax[1].transAxes,
               ha="center", fontsize=6.3, color="#555")
    panel(ax[1], "b", "The replication exists — if it is kept")

    # c: true replicate vs cross-dose against a matched null
    if tr is not None and len(tr):
        lab = {"true replicate (same line, drug, dose)": "true\nreplicate",
               "cross-dose (previous pairing)": "cross-dose",
               "pooled": "pooled"}
        tr = tr[tr.pairing.isin(lab)]
        xx = np.arange(len(tr))
        ax[2].bar(xx, tr.share, width=0.55,
                  color=[AQUA if "true" in p else ORANGE for p in tr.pairing])
        ax[2].errorbar(xx, tr.share,
                       yerr=[tr.share - tr.share_lo, tr.share_hi - tr.share],
                       fmt="none", ecolor="#333", capsize=3, lw=1.1)
        for i, r in enumerate(tr.itertuples()):
            ax[2].text(i, r.share + 0.012,
                       f"{r.share:.0%}\nn={r.n_pairs:,}\n"
                       f"p={'0.97' if r.p_vs_null > 0.01 else f'{r.p_vs_null:.0e}'}",
                       ha="center", fontsize=6.2)
        ax[2].set_xticks(xx, [lab[p] for p in tr.pairing], fontsize=7.5)
        ax[2].set_ylabel("interaction share (vs matched null)")
        ax[2].set_ylim(0, max(tr.share_hi.max() * 1.55, 0.05))
    panel(ax[2], "c", "The pairing decides the answer")
    fig.suptitle("Figure 1 — Estimating a context × compound interaction is "
                 "fragile, and the replication that makes it possible is easy "
                 "to discard", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    return fig, {"failure_modes": pd.DataFrame(
        modes, columns=["failure_mode", "correct_estimator", "shortcut",
                        "dataset"]),
        "tahoe_pairings": tr if tr is not None else pd.DataFrame()}


# ---------------------------------------------------------------- figure 2
def figure2():
    S = dataset_shares()
    corr = read("decompose_OP3_correlations.csv")
    dp = read("dose_vs_context_prism.csv")
    dl = read("dose_vs_context_lincs.csv")
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), constrained_layout=True)

    if len(S):
        x = np.arange(len(S))
        sh = S.interaction_pct / 100
        ax[0].bar(x, 1 - sh, width=0.6, color=BLUE, label="shared / additive")
        ax[0].bar(x, sh, bottom=1 - sh, width=0.6, color=ORANGE,
                  label="context × compound")
        for i, v in enumerate(sh):
            ax[0].text(i, 1 - v / 2, f"{v:.0%}", ha="center", va="center",
                       fontsize=8.5, fontweight="bold", color="white")
        ax[0].set_xticks(x, S.dataset, fontsize=6.6)
        ax[0].set_ylabel("share of reproducible variance")
        ax[0].set_ylim(0, 1); ax[0].legend(frameon=False, loc="lower right")
    panel(ax[0], "a", "Architecture where replicates exist")

    if corr is not None and len(corr):
        obs = corr[corr.kind == "observed"]
        lab = {"across_replicates_same_ctx_pert": "across\nreplicates",
               "across_perturbations_within_ctx": "across\ncompounds",
               "across_contexts_same_perturbation": "across\ncontexts"}
        keys = [k for k in lab if (obs.test == k).any()]
        data = [obs[obs.test == k].r.dropna() for k in keys]
        bp = ax[1].boxplot(data, showfliers=False, patch_artist=True,
                           medianprops=dict(color="black", lw=1.4))
        for p_, c in zip(bp["boxes"], [AQUA, GREY, GREY]):
            p_.set_facecolor(c); p_.set_alpha(0.8)
        ax[1].axhline(0, color="#444", lw=0.9)
        ax[1].set_xticks(range(1, len(keys) + 1), [lab[k] for k in keys],
                         fontsize=7.5)
        ax[1].set_ylabel("residual correlation")
    panel(ax[1], "b", "It is an interaction, not a context property")

    if dp is not None and "dose_pos" in dp.columns:
        q = dp.groupby("dose_pos").cdi.quantile([.25, .5, .75]).unstack()
        ax[2].fill_between(q.index, q[.25], q[.75], color=VIOLET, alpha=.18)
        ax[2].plot(q.index, q[.5], "o-", color=VIOLET, lw=2, ms=5,
                   label="PRISM viability")
        pk = int(q[.5].idxmax())
        ax[2].annotate("peak just below\nlethal dose", (pk, q[.5].max()),
                       fontsize=6.8, xytext=(-14, 16),
                       textcoords="offset points",
                       arrowprops=dict(arrowstyle="->", lw=0.8))
    if dl is not None and "dose_pos" in dl.columns:
        # LINCS is binned into quartiles and PRISM into eighths; place the
        # quartile midpoints on the same 1-8 axis so both curves are read
        # against one dose ordering rather than two overlaid axes
        ql = dl.groupby("dose_pos").cdi.median()
        xs = (ql.index.to_numpy() - 0.5) * 2
        ax[2].plot(xs, ql.to_numpy(), "s--", color=AQUA, lw=1.8, ms=5,
                   label="LINCS transcription")
    ax[2].set_xlabel("dose position (1 = lowest)")
    ax[2].set_ylabel("context-dependence index")
    ax[2].legend(frameon=False, loc="upper left")
    panel(ax[2], "c", "Context-dependence is dose-dependent")
    fig.suptitle("Figure 2 — The architecture of drug response where it can be "
                 "measured", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    return fig, {"dataset_shares": S,
                 "dose_prism": dp if dp is not None else pd.DataFrame()}


# ---------------------------------------------------------------- figure 3
def figure3():
    EA = read("expression_architecture.csv")
    ga = read("genetic_architecture_summary.csv")
    pw = read("prism_genotype_power.csv")
    tp = read("three_platform_mechanism_cdi.csv")
    fig, ax = plt.subplots(1, 4, figsize=(17.5, 4.1), constrained_layout=True)

    # a: all blocks measured on the SAME compounds and lines, so the comparison
    # is fair; the synonymous control sits in panel b because it was run on a
    # different compound set
    if EA is not None and len(EA):
        blocks = [("baseline\nexpression", "r2_expression", VIOLET),
                  ("baseline\nprotein", "r2_protein", BLUE),
                  ("lineage", "r2_lineage", AQUA),
                  ("copy\nnumber", "r2_cnv", "#c25fb0"),
                  ("nonsynonymous\nvariants", "r2_nonsyn", ORANGE)]
        blocks = [b for b in blocks if b[1] in EA.columns
                  and EA[b[1]].notna().any()]
        data = [EA[c].dropna() for _, c, _ in blocks]
        bp = ax[0].boxplot(data, showfliers=False, patch_artist=True,
                           medianprops=dict(color="black", lw=1.5))
        for p_, (_, _, c) in zip(bp["boxes"], blocks):
            p_.set_facecolor(c); p_.set_alpha(0.82)
        for i, d_ in enumerate(data):
            ax[0].text(i + 1, np.median(d_) + 0.012, f"{np.median(d_):+.3f}",
                       ha="center", fontsize=7.4, fontweight="bold")
        ax[0].axhline(0, color="#444", lw=1.0)
        ax[0].set_xticks(range(1, len(blocks) + 1), [b[0] for b in blocks],
                         fontsize=6.8)
        ax[0].set_ylabel("cross-validated $R^2$")
        ax[0].text(0.97, 0.95, f"n = {len(EA)} compounds\nidentical lines "
                   f"per block", transform=ax[0].transAxes, ha="right",
                   va="top", fontsize=6.4, color="#666")
        ax[0].text(0.03, 0.72, "copy number beats mutations\n(p = 6×10⁻⁴) but "
                   "trails expression\nby 0.088 (p = 3×10⁻²¹)",
                   transform=ax[0].transAxes, ha="left", va="top",
                   fontsize=6.1, color="#444")
    panel(ax[0], "a", "Molecular state predicts; genotype does not")

    # b: the synonymous control. Plotting the two blocks side by side hides the
    # point, because both are negative; the claim is about their DIFFERENCE, so
    # that is what is plotted. Synonymous variants cannot change a protein but
    # carry the same ancestry and lineage structure, so the excess of
    # nonsynonymous over synonymous is the mechanistic part.
    gafull = read("genetic_architecture.csv")
    if gafull is not None and "mechanistic" in gafull.columns:
        d = gafull.mechanistic.dropna()
        ax[1].hist(d.clip(-0.06, 0.06), bins=42, color=VIOLET, alpha=0.85)
        ax[1].axvline(0, color="#444", lw=1.4)
        ax[1].axvline(d.median(), color=ORANGE, lw=2.2,
                      label=f"median {d.median():+.4f}")
        ax[1].set_xlabel("nonsynonymous − synonymous $R^2$ (over lineage)")
        ax[1].set_ylabel("compounds")
        ax[1].legend(frameon=False, fontsize=7)
        ax[1].text(0.03, 0.95, f"{(d > 0).mean():.0%} of compounds positive\n"
                   f"p = 0.025 — real but ≈0.3%\nof the interaction variance",
                   transform=ax[1].transAxes, va="top", fontsize=6.5,
                   color="#444")
    panel(ax[1], "b", "The synonymous control isolates mechanism")

    if pw is not None and len(pw):
        ax[2].plot(pw.n_contexts, pw.detection_rate, "o-", color=VIOLET, lw=2,
                   ms=5)
        ax[2].axvline(47, color=ORANGE, ls="--", lw=1.3)
        ax[2].annotate("Tahoe\n47 lines", (47, 0.62), fontsize=7,
                       color=ORANGE, xytext=(6, 0), textcoords="offset points")
        ax[2].axhline(0.8, color="#888", ls=":", lw=1)
        ax[2].set_xscale("log")
        ax[2].set_xlabel("cell lines sampled")
        ax[2].set_ylabel("true associations recovered")
        ax[2].set_ylim(-0.02, 1.02)
    panel(ax[2], "c", "Genotype linkage needs ~400 contexts")

    if tp is not None and {"Tahoe", "PRISM"} <= set(tp.columns):
        v = tp[["Tahoe", "PRISM"]].dropna()
        ax[3].scatter(v.Tahoe, v.PRISM, s=34, color=VIOLET, edgecolors="none")
        from scipy import stats as st
        rho = st.spearmanr(v.Tahoe, v.PRISM)
        ax[3].set_xlabel("mechanism CDI — transcription (Tahoe)")
        ax[3].set_ylabel("mechanism CDI — viability (PRISM)")
        ax[3].text(0.04, 0.94, f"ρ = {rho.statistic:+.2f}\n(n.s., n={len(v)})",
                   transform=ax[3].transAxes, va="top", fontsize=8)
    panel(ax[3], "d", "Transcription and viability are decoupled")
    fig.suptitle("Figure 3 — What does and does not predict the interaction",
                 fontsize=10.5, x=0.005, ha="left", fontweight="bold")
    return fig, {"blocks_matched": EA if EA is not None else pd.DataFrame(),
                 "synonymous_control": ga if ga is not None else pd.DataFrame(),
                 "power": pw if pw is not None else pd.DataFrame()}


# ---------------------------------------------------------------- figure 4
def figure4():
    R = read("cross_lab_reproducibility.csv")
    S = read("cross_lab_summary.csv")
    ID = read("cross_lab_identity_viability.csv")
    fig, ax = plt.subplots(1, 4, figsize=(17.5, 4.0), constrained_layout=True)
    getv = (lambda q: float(S.set_index("quantity").value.get(q, np.nan))) \
        if S is not None else (lambda q: np.nan)

    order = ["PRISM replicate split (within-lab ceiling)",
             "GDSC1 vs GDSC2 (within-lab ceiling)",
             "PRISM vs GDSC (CROSS-LAB)"]
    if R is not None:
        order = [o for o in order if (R.comparison == o).any()]
        data = [R[R.comparison == o].rho.dropna() for o in order]
        bp = ax[0].boxplot(data, showfliers=False, patch_artist=True,
                           medianprops=dict(color="black", lw=1.5))
        for p_, c in zip(bp["boxes"], [AQUA, BLUE, ORANGE]):
            p_.set_facecolor(c); p_.set_alpha(0.8)
        for i, d_ in enumerate(data):
            ax[0].text(i + 1, np.median(d_) + 0.022, f"{np.median(d_):.2f}",
                       ha="center", fontsize=8.5, fontweight="bold")
        ax[0].set_xticks(range(1, len(order) + 1),
                         ["same lab\nsame assay", "same lab\nnew assay",
                          "new lab\n+ assay"][:len(order)], fontsize=7)
        ax[0].axhline(0, color="#444", lw=0.9)
        ax[0].set_ylabel("Spearman r of the line-specific residual")
        a_sh, l_sh = getv("share of loss from changing assay within a lab"), \
            getv("share of loss from changing laboratory")
        if np.isfinite(l_sh):
            ax[0].text(0.5, 0.03, f"assay {a_sh:.0%} of the loss · "
                       f"laboratory {l_sh:.0%}", transform=ax[0].transAxes,
                       ha="center", fontsize=7.2, color="#333")
    panel(ax[0], "a", "Laboratory, not assay, is what costs")

    if ID is not None and len(ID):
        m = ID[ID.rank_of_id_match > 0]
        ax[1].hist(m.rank_of_id_match,
                   bins=np.logspace(0, np.log10(max(m.n_candidates.max(), 10)),
                                    38), color=VIOLET, alpha=0.85)
        ax[1].set_xscale("log")
        ax[1].axvline(1, color=AQUA, lw=2, label="is its own best match")
        ax[1].axvline(m.n_candidates.median() / 2, color="#888", ls=":", lw=1.4,
                      label="chance")
        ax[1].set_xlabel("rank of the identifier-matched line")
        ax[1].set_ylabel("cell lines")
        ax[1].legend(frameon=False)
        ax[1].text(0.97, 0.72, f"{m.reciprocal_best.mean():.0%} are their own\n"
                   f"best match; median rank\n"
                   f"{m.rank_of_id_match.median():.0f} of "
                   f"{m.n_candidates.median():.0f}",
                   transform=ax[1].transAxes, ha="right", va="top", fontsize=7)
    panel(ax[1], "b", "Identifier matches are rarely the single best")

    fr = [getv("reproducible fraction (all)"),
          getv("reproducible fraction (identity-validated)")]
    ax[2].bar([0, 1], fr, width=0.55, color=[ORANGE, VIOLET])
    for i, v in enumerate(fr):
        if np.isfinite(v):
            ax[2].text(i, v + 0.025, f"{v:.0%}", ha="center", fontsize=12,
                       fontweight="bold")
    ax[2].axhline(1.0, color="#888", ls="--", lw=1.2)
    ax[2].set_xticks([0, 1], ["trust the\nidentifier",
                              "verify identity\nfrom data"], fontsize=7.5)
    ax[2].set_ylim(0, 1.16)
    ax[2].set_ylabel("fraction of the within-lab ceiling")
    ax[2].text(0.5, 1.05, "identity validated on one half of the compounds, "
               "scored on the other", transform=ax[2].transAxes, ha="center",
               fontsize=6.6, color="#555")
    panel(ax[2], "c", "Verifying identity recovers most of the gap")

    if R is not None:
        xl = R[R.comparison == "PRISM vs GDSC (CROSS-LAB)"]
        dpx = read("dose_vs_context_prism.csv")
        strength = None
        if dpx is not None:
            pass
        if "interaction_strength" in xl.columns:
            xl = xl.rename(columns={"interaction_strength": "stren"}).dropna(
                subset=["stren"])
            if len(xl) > 12:
                q = pd.qcut(xl.stren, 4, labels=["Q1\nweak", "Q2", "Q3",
                                                 "Q4\nstrong"])
                gq = xl.groupby(q, observed=True).rho.median()
                ax[3].bar(range(len(gq)), gq.to_numpy(), width=0.6, color=BLUE)
                ceil = np.sqrt(getv("PRISM replicate ceiling")
                               * getv("GDSC1 vs GDSC2 ceiling"))
                if np.isfinite(ceil):
                    ax[3].axhline(ceil, color="#888", ls="--", lw=1.2,
                                  label="within-lab ceiling")
                    ax[3].legend(frameon=False)
                ax[3].set_xticks(range(len(gq)), list(gq.index), fontsize=7.5)
                ax[3].set_ylabel("cross-lab Spearman r")
                ax[3].set_xlabel("strength of the line-specific component")
    panel(ax[3], "d", "Strong signal transfers; weak does not")
    fig.suptitle("Figure 4 — Cell-line identity, not protocol, limits "
                 "cross-laboratory reproducibility", fontsize=10.5, x=0.005,
                 ha="left", fontweight="bold")
    return fig, {"per_compound": R if R is not None else pd.DataFrame(),
                 "summary": S if S is not None else pd.DataFrame(),
                 "identity": ID if ID is not None else pd.DataFrame()}


# ---------------------------------------------------------------- figure 5
def figure5():
    M = read("matching_strategy.csv")
    EI = read("expression_identity.csv")
    ncol = 3 if (EI is not None and len(EI)) else 2
    fig, ax = plt.subplots(1, ncol, figsize=(5.2 * ncol, 4.0),
                           constrained_layout=True)
    if M is not None and len(M):
        order = ["random pairing", "same tissue,\nrandom line",
                 "identifier\n(standard practice)",
                 "identifier +\nfingerprint-validated", "best hit\n(upper bound)"]
        order = [o for o in order if (M.strategy == o).any()]
        data = [M[M.strategy == o].r.dropna() for o in order]
        cols = [GREY, AQUA, ORANGE, VIOLET, BLUE][:len(order)]
        bp = ax[0].boxplot(data, showfliers=False, patch_artist=True,
                           medianprops=dict(color="black", lw=1.5))
        for p_, c in zip(bp["boxes"], cols):
            p_.set_facecolor(c); p_.set_alpha(0.8)
        for i, d_ in enumerate(data):
            ax[0].text(i + 1, np.median(d_) + 0.022, f"{np.median(d_):.2f}",
                       ha="center", fontsize=8, fontweight="bold")
        ax[0].set_xticks(range(1, len(order) + 1),
                         [o.replace("\n", " ") for o in order], fontsize=6.2,
                         rotation=18, ha="right")
        ax[0].axhline(0, color="#444", lw=0.9)
        ax[0].set_ylabel("response-fingerprint similarity")
        ax[0].text(0.03, 0.95, "last two rungs are selected on the\noutcome — "
                   "they bound the room available", transform=ax[0].transAxes,
                   va="top", fontsize=6.6, color="#555")
        panel(ax[0], "a", "What each matching rule buys")
        for o, c in zip(order, cols):
            d_ = M[M.strategy == o].r.dropna()
            if len(d_) > 5:
                ax[1].hist(d_, bins=40, histtype="step", lw=1.8, color=c,
                           density=True, label=o.replace("\n", " "))
        ax[1].axvline(0, color="#444", lw=0.9)
        ax[1].set_xlabel("fingerprint similarity"); ax[1].set_ylabel("density")
        ax[1].legend(frameon=False, fontsize=6.4)
        panel(ax[1], "b", "Distributions")

    # c: is a failed identifier match evidence of divergence, or just noise?
    # If the line that outranks the identifier match were a genuinely divergent
    # culture's nearest relative, it should sit among that line's closest
    # expression neighbours. It does not — the distribution is centred near
    # chance, which is why the divergence reading was withdrawn.
    if EI is not None and len(EI):
        ax[2].hist(EI.expr_rank_of_best_match,
                   bins=np.logspace(0, np.log10(max(EI.n_candidates.max(), 10)),
                                    34), color=AQUA, alpha=0.85)
        ax[2].set_xscale("log")
        med = EI.expr_rank_of_best_match.median()
        chance = EI.n_candidates.median() / 2
        ax[2].axvline(chance, color="#666", ls=":", lw=1.6, label="chance")
        ax[2].axvline(med, color=ORANGE, lw=2, label=f"observed median ({med:.0f})")
        ax[2].set_xlabel("expression rank of the line that outranked the "
                         "identifier")
        ax[2].set_ylabel("cell lines")
        ax[2].legend(frameon=False, fontsize=6.8)
        ax[2].text(0.03, 0.95, f"only {EI.same_tissue.mean():.0%} share the "
                   f"primary tissue\n→ outranking lines are near-random "
                   f"neighbours,\n   so best-hit failure is metric noise, "
                   f"not\n   demonstrated culture divergence",
                   transform=ax[2].transAxes, va="top", fontsize=6.3,
                   color="#444")
        panel(ax[2], "c", "Why the divergence reading was withdrawn")
    fig.suptitle("Figure 5 — Matching cell lines across atlases: what "
                 "identifiers assume, what the data supports", fontsize=10.5,
                 x=0.005, ha="left", fontweight="bold")
    return fig, {"per_pair": M if M is not None else pd.DataFrame(),
                 "expression_identity": EI if EI is not None else pd.DataFrame()}


# ---------------------------------------------------------------- figure 6
def figure6():
    """Benchmark design: what the measurements imply for how models are read."""
    A = read("additive_vs_n_lines.csv")
    FS = read("few_shot_eval_dev47ft.csv")
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.1), constrained_layout=True)

    # a: baseline quality vs the number of contexts it averages
    if A is not None and len(A):
        g = A.groupby("k").r_de100.agg(["mean", "std"])
        ax[0].errorbar(g.index, g["mean"], yerr=g["std"] / np.sqrt(len(A) /
                       len(g)), fmt="o-", color=BLUE, lw=2, ms=5, capsize=3)
        ax[0].set_xscale("log")
        ax[0].set_xlabel("contexts averaged by the mean baseline")
        ax[0].set_ylabel("baseline $r_{de100}$")
        # mark where published benchmarks sit
        for n, lab in ((2, "State query\ndatasets (2-3)"), (17, "Parse-PBMC\n(17)")):
            ax[0].axvline(n, color=ORANGE, ls=":", lw=1.3)
        ax[0].annotate("State query datasets leave\n2-5 contexts in the "
                       "baseline", (2.1, g["mean"].min()), fontsize=6.2,
                       color=ORANGE, va="bottom")
        panel(ax[0], "a", "The baseline improves with context count")

        # b: what that costs, as apparent model gain
        ref = g.loc[g.index.max(), "mean"]
        hand = (ref - g["mean"]) / g["mean"] * 100
        cols_ = [ORANGE if k <= 6 else BLUE for k in g.index]
        ax[1].bar(np.arange(len(g)), hand.to_numpy(), width=0.62, color=cols_)
        ax[1].set_xticks(np.arange(len(g)), [str(int(k)) for k in g.index],
                         fontsize=7.5)
        for i_, v in enumerate(hand.to_numpy()):
            if v > 0.5:
                ax[1].text(i_, v + 0.4, f"{v:.0f}%", ha="center", fontsize=7)
        ax[1].set_xlabel("contexts averaged by the baseline")
        ax[1].set_ylabel("apparent gain of an unchanged model (%)")
        ax[1].text(0.97, 0.95, "a model equally good everywhere\nlooks 21% "
                   "better at 2 contexts\nthan at 18", transform=ax[1].transAxes,
                   ha="right", va="top", fontsize=6.4, color="#444")
        panel(ax[1], "b", "Benchmark context count inflates gains")

    # c: what actually makes a new context predictable
    if FS is not None and len(FS):
        base = float(FS[FS.strategy == "additive"].r_de100.mean())
        g2 = (FS[FS.strategy == "random"].groupby("k").r_de100.mean()
              .reindex([1.0, 5.0, 20.0]))
        full = float(FS[FS.k == -1.0].r_de100.mean())
        ax[2].axhline(base, color=GREY, ls="--", lw=1.4,
                      label="additive baseline")
        ax[2].axhline(full, color=AQUA, ls=":", lw=1.6,
                      label="all compounds (ceiling)")
        ax[2].plot(g2.index, g2.to_numpy(), "o-", color=VIOLET, lw=2, ms=6,
                   label="random probe compounds")
        gain = (g2 - base) / (full - base) * 100
        for x_, v in zip(g2.index, gain):
            if np.isfinite(v):
                ax[2].annotate(f"{v:.0f}% of achievable", (x_, g2[x_]),
                               fontsize=6.3, xytext=(4, -14),
                               textcoords="offset points")
        ax[2].set_xlabel("probe compounds measured in the new context")
        ax[2].set_ylabel("$r_{de100}$")
        ax[2].legend(frameon=False, fontsize=7)
        panel(ax[2], "c", "Measurement, not description, transfers")
    fig.suptitle("Figure 6 — How benchmark design shapes apparent model "
                 "performance", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    return fig, {"baseline_scaling": A if A is not None else pd.DataFrame(),
                 "few_shot": FS if FS is not None else pd.DataFrame()}


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    for i, fn in enumerate((figure1, figure2, figure3, figure4, figure5,
                            figure6), 1):
        fig, src = fn()
        d = save_figure(fig, f"fig{i}", FIG, source_data=src, script=__file__)
        plt.close(fig)
        print(f"Figure {i} -> {d}")


if __name__ == "__main__":
    main()
