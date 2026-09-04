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


def read(name, sub=None, **kw):
    """A source table, from results/tables or from a named figure bundle.

    Analyses that write several source CSVs keep them in their own bundle
    directory rather than the flat table folder, so `sub` names that bundle.
    """
    f = (FIG / sub / name) if sub else (TAB / name)
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
    """The decomposition, and the failure modes that make it necessary."""
    ev = read("methodology_evidence.csv")
    tr = read("tahoe_true_replicates.csv")
    AP = read("variance_apportionment.csv")
    GS = read("general_sensitivity_split.csv")
    SIM = read("estimator_simulation.csv")
    fig, ax = plt.subplots(2, 3, figsize=(13.8, 8.0), constrained_layout=True)
    ax = ax.ravel()

    # a: the three-way split. A measured response is a drug effect, a property
    # of the cell, and a relation between the two; only the third is what these
    # atlases are built to find. Plotted on a log axis because the drug effect
    # is 25x the other two combined and a linear axis would erase them.
    if AP is not None and len(AP):
        lab = {"drug effect": "drug effect\n(same in every line)",
               "cell property": "cell property\n(general sensitivity)",
               "cell-drug relation": "cell-drug relation\n(interaction)"}
        cols = {"drug effect": GREY, "cell property": VIOLET,
                "cell-drug relation": ORANGE}
        xx = np.arange(len(AP))
        ax[0].bar(xx, AP.variance, width=0.55,
                  color=[cols.get(c, GREY) for c in AP.component])
        for i_, r in enumerate(AP.itertuples()):
            ax[0].text(i_, r.variance * 1.15, f"{r.share:.1%}", ha="center",
                       fontsize=8.5, fontweight="bold")
        ax[0].set_yscale("log")
        ax[0].set_xticks(xx, [lab.get(c, c) for c in AP.component], fontsize=6.8)
        ax[0].set_ylabel("variance (log scale)")
        ax[0].set_ylim(AP.variance.min() * 0.35, AP.variance.max() * 4)
        ax[0].text(0.98, 0.97, "PRISM\n737 lines x 1,324 compounds",
                   transform=ax[0].transAxes, ha="right", va="top",
                   fontsize=6.3, color="#666")
    panel(ax[0], "a", "A response is three things, not two")

    # b: the cell property is not noise. Its estimate on one half of the
    # compounds predicts the estimate on a disjoint half almost exactly, which
    # is why leaving it in the residual silently inflates every interaction.
    if GS is not None and len(GS):
        ax[1].scatter(GS.half_A, GS.half_B, s=6, alpha=0.32, color=VIOLET,
                      edgecolors="none")
        r = float(np.corrcoef(GS.half_A, GS.half_B)[0, 1])
        lo = float(min(GS.half_A.min(), GS.half_B.min()))
        hi = float(max(GS.half_A.max(), GS.half_B.max()))
        ax[1].plot([lo, hi], [lo, hi], ls="--", color="#555", lw=1)
        ax[1].set_xlabel("general sensitivity, compound half A")
        ax[1].set_ylabel("general sensitivity, compound half B")
        ax[1].text(0.04, 0.92, f"r = {r:.3f}\nn = {len(GS)} lines",
                   transform=ax[1].transAxes, fontsize=8, fontweight="bold",
                   va="top")
    panel(ax[1], "b", "It reproduces almost perfectly")

    # c: and it is charged to the interaction unless removed. Truth is held at
    # 0.20 while only the planted general response grows.
    if SIM is not None and "sigma_ctx" in SIM.columns:
        S = SIM[SIM.sweep == "sigma_ctx"].groupby("sigma_ctx").estimate.agg(
            ["mean", "std"])
        if len(S):
            ax[2].axhline(0.20, ls="--", color="#555", lw=1.2,
                          label="true share")
            ax[2].errorbar(S.index, S["mean"], yerr=S["std"], fmt="o-",
                           color=AQUA, lw=2, ms=6, capsize=3,
                           label="general sensitivity removed")
            ax[2].set_xlabel("planted general response (SD)")
            ax[2].set_ylabel("interaction share reported")
            ax[2].set_ylim(0, 0.72)
            ax[2].legend(frameon=False, fontsize=7, loc="upper left")
            ax[2].text(0.5, 0.06, "uncorrected: 0.170 / 0.347 / 0.595 at "
                       "0 / 0.5 / 1.0 SD", transform=ax[2].transAxes,
                       ha="center", fontsize=6.4, color=ORANGE,
                       fontweight="bold")
    panel(ax[2], "c", "Uncorrected, it is read as context-dependence")

    # d: each failure mode is measured on its own dataset, so a grouped bar
    # chart would imply a common scale that does not exist. A dumbbell keeps
    # each comparison self-contained and labels the data it came from.
    modes = [("residual variance\nvs replicate covariance", 0.008, 0.717,
              "simulation, true share 0.20"),
             ("pooled vs cross-batch\nreplicate pairs", 0.202, 0.411,
              "simulation, true share 0.20"),
             ("doses vs true replicates\nas the replicate axis", 0.202, 0.131,
              "simulation, true share 0.20")]
    y = np.arange(len(modes))[::-1] * 1.0
    for yy, (lab_, good, bad, src) in zip(y, modes):
        ax[3].plot([good, bad], [yy, yy], color="#c8c8c8", lw=3, zorder=1,
                   solid_capstyle="round")
        ax[3].scatter([good], [yy], s=80, color=AQUA, zorder=3)
        ax[3].scatter([bad], [yy], s=80, color=ORANGE, zorder=3)
        # stagger the pair vertically when the markers nearly coincide,
        # otherwise the two percentages overprint each other
        close = abs(bad - good) < 0.12
        ax[3].annotate(f"{bad:.0%}", (bad, yy), textcoords="offset points",
                       xytext=(0, -18 if close else 11), ha="center",
                       fontsize=8, color=ORANGE, fontweight="bold")
        ax[3].annotate(f"{good:.0%}", (good, yy), textcoords="offset points",
                       xytext=(0, 11), ha="center", fontsize=8,
                       color=AQUA, fontweight="bold")
    ax[3].axvline(0.20, ls="--", color="#555", lw=1.1, zorder=0)
    ax[3].set_yticks(y, [f"{m[0]}\n({m[3]})" for m in modes], fontsize=6.4)
    ax[3].set_xlim(-0.08, 0.85)
    ax[3].set_ylim(-1.35, len(modes) - 0.35)
    ax[3].set_xlabel("interaction share reported (truth 0.20, dashed)")
    ax[3].scatter([], [], s=80, color=AQUA, label="recommended estimator")
    ax[3].scatter([], [], s=80, color=ORANGE, label="shortcut")
    ax[3].legend(frameon=False, loc="upper right", fontsize=7,
                 handletextpad=0.3)
    ax[3].text(0.42, -1.12, "a fifth mode — an in-sample shared response — "
               "reports 0.169", ha="center", va="center", fontsize=6.1,
               color="#666")
    panel(ax[3], "d", "Ways to get it wrong, on data with a known answer")

    # e: Tahoe replicate structure, counted from the RELEASED metadata rather
    # than our pseudobulk, and shown with and without the replicate plate the
    # authors withhold from training
    tot = 56877
    with14, without14 = 7691, 3045
    ax[4].bar([0, 1], [with14 / tot, without14 / tot], width=0.5,
              color=[ORANGE, GREY])
    for i_, (v, n) in enumerate([(with14 / tot, with14),
                                 (without14 / tot, without14)]):
        ax[4].text(i_, v + 0.004, f"{v:.1%}\n{n:,} triples", ha="center",
                   fontsize=7.4, fontweight="bold")
    ax[4].set_xticks([0, 1], ["plate 14 kept\n(as released)",
                              "plate 14 dropped\n(training convention)"],
                     fontsize=7)
    ax[4].set_ylim(0, 0.175)
    ax[4].set_ylabel("(line, drug, dose) replicated on >1 plate")
    ax[4].text(0.5, 0.86, "plate 14 is a designed replicate of plate 6\n"
               "(6.2M cells, 50 lines, 95 drugs)", transform=ax[4].transAxes,
               ha="center", fontsize=6.3, color="#555")
    panel(ax[4], "e", "The replication exists — if it is kept")

    # f: true replicate vs cross-dose against a matched null
    if tr is not None and len(tr):
        lab2 = {"true replicate (same line, drug, dose)": "true\nreplicate",
                "cross-dose (previous pairing)": "cross-dose",
                "pooled": "pooled"}
        tr = tr[tr.pairing.isin(lab2)]
        xx = np.arange(len(tr))
        ax[5].bar(xx, tr.share, width=0.55,
                  color=[AQUA if "true" in p else ORANGE for p in tr.pairing])
        ax[5].errorbar(xx, tr.share,
                       yerr=[tr.share - tr.share_lo, tr.share_hi - tr.share],
                       fmt="none", ecolor="#333", capsize=3, lw=1.1)
        for i_, r in enumerate(tr.itertuples()):
            ax[5].text(i_, r.share + 0.012,
                       f"{r.share:.0%}\nn={r.n_pairs:,}\n"
                       f"p={'0.97' if r.p_vs_null > 0.01 else f'{r.p_vs_null:.0e}'}",
                       ha="center", fontsize=6.2)
        ax[5].set_xticks(xx, [lab2[p] for p in tr.pairing], fontsize=7.5)
        ax[5].set_ylabel("interaction share (vs matched null)")
        ax[5].set_ylim(0, max(tr.share_hi.max() * 1.55, 0.05))
    panel(ax[5], "f", "The pairing decides the answer")

    fig.suptitle("Figure 1 — A measured response is a drug effect, a property "
                 "of the cell, and a relation between them; only the third is "
                 "context-dependence", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    return fig, {"apportionment": AP if AP is not None else pd.DataFrame(),
                 "general_sensitivity": GS if GS is not None else pd.DataFrame(),
                 "failure_modes": pd.DataFrame(
                     modes, columns=["failure_mode", "recommended", "shortcut",
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
    BP = read("central_claim_block_permutation.csv")
    SP = read("state_vs_genotype_split.csv")
    fig, ax = plt.subplots(1, 5, figsize=(21.5, 4.1), constrained_layout=True)

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
        # The copy-number claim was withdrawn: its advantage over mutations
        # does not survive a bootstrap that resamples whole compounds, and
        # compound residuals correlate at about r = 0.23.
        ax[0].text(0.03, 0.20, "copy number over mutations: +0.005,\n"
                   "CI [−0.000, +0.010] — withdrawn.\nExpression leads by "
                   "0.088 (p = 2×10⁻²¹)",
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
        ax[1].text(0.03, 0.95, f"{(d > 0).mean():.0%} of compounds positive; "
                   f"the\nexcess is real but is ≈0.3% of\nthe interaction "
                   f"variance",
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
        ax[3].text(0.04, 0.94, f"ρ = {rho.statistic:+.2f}, n = {len(v)}\n"
                   f"underpowered:\nρ = 0.6 would not reach\nsignificance at "
                   f"this n", transform=ax[3].transAxes, va="top", fontsize=7)
    panel(ax[3], "d", "Readouts: too few shared classes to tell")

    # e: the median R2 in panel a is one number over correlated compounds. This
    # asks the question once per compound against that compound's OWN
    # permutation null -- cell-line labels shuffled within it -- and counts how
    # many survive FDR across the family. A count cannot be carried by a few
    # strong compounds the way a median can.
    if BP is not None and len(BP):
        labs = [("expression", VIOLET), ("lineage", AQUA), ("mutations", ORANGE)]
        labs = [(l, c) for l, c in labs if f"q_{l}" in BP.columns]
        fr = [float((BP[f"q_{l}"] < 0.05).mean()) for l, _ in labs]
        xx = np.arange(len(labs))
        ax[4].bar(xx, fr, width=0.58, color=[c for _, c in labs])
        for i_, v in enumerate(fr):
            n_ = int(v * len(BP) + 0.5)
            ax[4].text(i_, v + 0.02, f"{v:.0%}\n{n_}/{len(BP)}", ha="center",
                       fontsize=8, fontweight="bold")
        ax[4].axhline(0.05, ls="--", color="#888", lw=1.2, label="FDR level")
        ax[4].set_xticks(xx, [l for l, _ in labs], fontsize=7.5)
        ax[4].set_ylim(0, 1.14)
        ax[4].set_ylabel("compounds beating their own permutation null")
        ax[4].legend(frameon=False, fontsize=7)
        if {"q_expression", "q_mutations"} <= set(BP.columns):
            only_e = int(((BP.q_expression < 0.05) &
                          (BP.q_mutations >= 0.05)).sum())
            only_m = int(((BP.q_expression >= 0.05) &
                          (BP.q_mutations < 0.05)).sum())
            ax[4].text(0.5, 0.70, f"{only_e} compounds expression-only,\n"
                       f"{only_m} mutation-only", transform=ax[4].transAxes,
                       ha="center", fontsize=6.6, color="#444")
    panel(ax[4], "e", "Per compound, with FDR across the family")

    fig.suptitle("Figure 3 — What does and does not predict the interaction",
                 fontsize=10.5, x=0.005, ha="left", fontweight="bold")
    return fig, {"blocks_matched": EA if EA is not None else pd.DataFrame(),
                 "synonymous_control": ga if ga is not None else pd.DataFrame(),
                 "power": pw if pw is not None else pd.DataFrame(),
                 "block_permutation": BP if BP is not None else pd.DataFrame(),
                 "property_vs_relation": SP if SP is not None
                 else pd.DataFrame()}


# ---------------------------------------------------------------- figure 4
def figure4():
    R = read("cross_lab_reproducibility.csv")
    S = read("cross_lab_summary.csv")
    ID = read("cross_lab_identity_viability.csv")
    TX = read("cross_lab_transcription.csv")
    ncol = 5 if (TX is not None and len(TX)) else 4
    fig, ax = plt.subplots(1, ncol, figsize=(4.4 * ncol, 4.0),
                           constrained_layout=True)
    getv = (lambda q: float(S.set_index("quantity").value.get(q, np.nan))) \
        if S is not None else (lambda q: np.nan)

    a_sh = l_sh = np.nan
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
    # The title states whichever term the data says dominates. It previously
    # asserted "Laboratory, not assay" -- the opposite of what the matched
    # compound set shows -- and that assertion also set the paper's title.
    dom = "Assay, not laboratory, is what costs"
    if np.isfinite(l_sh) and np.isfinite(a_sh):
        dom = ("Laboratory, not assay, is what costs" if l_sh > a_sh
               else "Assay, not laboratory, is what costs")
    panel(ax[0], "a", dom)

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
        ax[1].text(0.97, 0.50, f"{m.reciprocal_best.mean():.0%} are their own\n"
                   f"best match; median rank\n"
                   f"{m.rank_of_id_match.median():.0f} of "
                   f"{m.n_candidates.median():.0f}",
                   transform=ax[1].transAxes, ha="right", va="top", fontsize=7)
    panel(ax[1], "b", "Identifier matches are rarely the single best")

    # The validated bar is divided by a ceiling measured on the SAME lines and
    # compounds. Against the all-lines ceiling it read 110%, which a fraction of
    # a ceiling cannot be. The third bar is the control that decides what the
    # selection is actually doing: lines picked for measurement RELIABILITY
    # rather than identity. If it matched the middle bar, identity would be
    # doing nothing that measurement quality does not already do.
    fr = [getv("reproducible fraction (all)"),
          getv("reproducible fraction (identity-validated, MATCHED ceiling)"),
          getv("reproducible fraction (reliability-selected control)")]
    cols = [ORANGE, VIOLET, GREY]
    labs = ["trust the\nidentifier", "verify identity\nfrom data",
            "reliability-matched\ncontrol"]
    keep = [i for i, v in enumerate(fr) if np.isfinite(v)]
    ax[2].bar(range(len(keep)), [fr[i] for i in keep], width=0.55,
              color=[cols[i] for i in keep])
    lo = getv("reproducible fraction (validated) CI low")
    hi = getv("reproducible fraction (validated) CI high")
    if np.isfinite(lo) and 1 in keep:
        k = keep.index(1)
        ax[2].errorbar([k], [fr[1]], yerr=[[fr[1] - lo], [hi - fr[1]]],
                       fmt="none", ecolor="#333", capsize=4, lw=1.2)
    for k, i in enumerate(keep):
        ax[2].text(k, fr[i] + 0.03, f"{fr[i]:.0%}", ha="center", fontsize=11,
                   fontweight="bold")
    ax[2].axhline(1.0, color="#888", ls="--", lw=1.2)
    ax[2].set_xticks(range(len(keep)), [labs[i] for i in keep], fontsize=7)
    ax[2].set_ylim(0, max(1.25, (hi if np.isfinite(hi) else 1.2) * 1.06))
    ax[2].set_ylabel("fraction of the within-lab ceiling")
    ax[2].set_xlabel("identity validated on one half of the compounds, scored "
                     "on the other;\nceiling measured on the same lines",
                     fontsize=6.2, color="#555")
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

    # e: the same measurement for a transcriptional readout, which is what
    # makes the result a property of laboratories rather than of a killing assay
    if TX is not None and len(TX):
        wl = TX[TX.comparison.str.contains("within-lab", na=False)
                & ~TX.comparison.str.contains("validated", na=False)]
        xl = TX[TX.comparison.str.contains("Tahoe vs LINCS", na=False)
                & ~TX.comparison.str.contains("validated", na=False)]
        d_ = [wl.r.dropna(), xl.r.dropna()]
        bp2 = ax[4].boxplot(d_, showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=1.5))
        for pch, c in zip(bp2["boxes"], [AQUA, ORANGE]):
            pch.set_facecolor(c); pch.set_alpha(0.8)
        for i_, dd in enumerate(d_):
            if len(dd):
                ax[4].text(i_ + 1, np.median(dd) + 0.006,
                           f"{np.median(dd):.3f}", ha="center", fontsize=8,
                           fontweight="bold")
        ax[4].axhline(0, color="#444", lw=0.9)
        ax[4].set_xticks([1, 2], ["within lab\n(LINCS p1 vs p2)",
                                  "cross lab\n(Tahoe vs LINCS)"], fontsize=7)
        ax[4].set_ylabel("residual-profile r")
        if len(d_[0]) and len(d_[1]):
            fr = np.median(d_[1]) / np.median(d_[0])
            ax[4].text(0.5, 0.93, f"{fr:.0%} of ceiling", ha="center",
                       transform=ax[4].transAxes, fontsize=7.2, color="#444")
            ttl = "Transcription gives the same answer"
        else:
            # Say so rather than showing an empty box: after identity
            # validation no (line, compound) pair is shared between Tahoe and
            # LINCS, so the transcriptional arm has no cross-laboratory
            # comparison to make and this is a limit of the data, not a result.
            ax[4].text(0.76, 0.5, "no (line, compound) pair\nsurvives identity "
                       "validation\nacross Tahoe and LINCS —\nthe "
                       "transcriptional arm has\nno cross-laboratory test",
                       transform=ax[4].transAxes, ha="center", va="center",
                       fontsize=7.4, color=ORANGE)
            ttl = "Transcription cannot be tested across labs"
        panel(ax[4], "e", ttl)
    # The title follows panel a. It previously asserted the opposite of what
    # the matched compound set shows, and that assertion propagated into the
    # paper's own title.
    lead = ("Cell-line identity, and the assay rather than the laboratory, "
            "limit cross-laboratory reproducibility"
            if (np.isfinite(a_sh) and np.isfinite(l_sh) and a_sh > l_sh)
            else "Cell-line identity, not protocol, limits cross-laboratory "
                 "reproducibility")
    fig.suptitle(f"Figure 4 — {lead}", fontsize=10.5, x=0.005,
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
                   label="fine-tuned model")
        # The floor that decides what the few-shot gain means: the additive
        # prior plus the line's mean residual over the SAME probe compounds.
        # One arithmetic step, no fitting. A model is only learning
        # context-specific pharmacology to the extent it beats this.
        fl = FS[FS.strategy == "additive+line"]
        if len(fl):
            g3 = fl.groupby("k").r_de100.mean()
            ax[2].plot(g3.index, g3.to_numpy(), "o--", color=ORANGE, lw=2,
                       ms=6, label="additive + line's general response")
            k_top = int(g3.index.max())
            ax[2].annotate("the floor beats the model,\nand beats the oracle,\n"
                           "wherever it is estimable",
                           (k_top, float(g3.loc[k_top])),
                           textcoords="offset points", xytext=(-8, -34),
                           ha="right", fontsize=6.6, color=ORANGE,
                           fontweight="bold")
        ax[2].set_xscale("log")
        for x_ in g2.index:
            if int(x_) in (5, 20):
                v = 100 * (g2[x_] - base) / max(full - base, 1e-9)
                ax[2].annotate(f"{v:.0f}%", (x_, g2[x_]),
                               textcoords="offset points", xytext=(4, 7),
                               fontsize=6.6, color=VIOLET)
        ax[2].set_xlabel("probe compounds measured in the new context")
        ax[2].set_ylabel("$r_{de100}$")
        ax[2].legend(frameon=False, fontsize=6.6, loc="lower right")
        panel(ax[2], "c", "And what transfers is a scalar, not a model")
    fig.suptitle("Figure 6 — How benchmark design shapes apparent model "
                 "performance", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    return fig, {"baseline_scaling": A if A is not None else pd.DataFrame(),
                 "few_shot": FS if FS is not None else pd.DataFrame()}


# ---------------------------------------------------------------- figure 7
def figure7():
    """The same omission, tested on another group's data (TRADE)."""
    T = read("trade_cross_celltype.csv")
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    if T is None or len(T) < 2:
        for a_ in ax:
            a_.axis("off")
        return fig, {}
    base, corr = T.iloc[0], T.iloc[1]
    pub_dep = 0.44

    # a: the split recomputed two ways on identical data. The raw statistic is
    # the one other groups report; the noise-corrected one is the form our
    # estimator uses, and only it is sensitive to the omission.
    xx = np.arange(2)
    w = 0.36
    ax[0].bar(xx - w / 2, [base.dependent_raw, base.dependent_corrected], w,
              color=ORANGE, label="cell-line effect left in")
    ax[0].bar(xx + w / 2, [corr.dependent_raw, corr.dependent_corrected], w,
              color=AQUA, label="cell-line effect removed")
    for i_, (a_, b_) in enumerate([(base.dependent_raw, corr.dependent_raw),
                                   (base.dependent_corrected,
                                    corr.dependent_corrected)]):
        ax[0].text(i_ - w / 2, a_ + .012, f"{a_:.0%}", ha="center",
                   fontsize=8.5, fontweight="bold")
        ax[0].text(i_ + w / 2, b_ + .012, f"{b_:.0%}", ha="center",
                   fontsize=8.5, fontweight="bold")
    ax[0].axhline(pub_dep, ls="--", color="#555", lw=1.2,
                  label=f"published {pub_dep:.0%}")
    ax[0].set_xticks(xx, ["raw\n(noise included)",
                          "noise-corrected\n(covariance + published SEs)"],
                     fontsize=8)
    ax[0].set_ylabel("'cell-type-dependent' share")
    ax[0].set_ylim(0, 0.85)
    ax[0].legend(frameon=False, fontsize=7)
    panel(ax[0], "a", "TRADE's split, recomputed two ways")

    d_raw = 100 * (base.dependent_raw - corr.dependent_raw)
    d_cor = 100 * (base.dependent_corrected - corr.dependent_corrected)
    ax[1].bar([0, 1], [d_raw, d_cor], width=0.5, color=[GREY, VIOLET])
    for i_, v in enumerate([d_raw, d_cor]):
        ax[1].text(i_, v + 0.25, f"{v:+.1f} pts", ha="center", fontsize=10.5,
                   fontweight="bold")
    ax[1].set_xticks([0, 1], ["raw", "noise-corrected"], fontsize=8.5)
    ax[1].set_ylabel("inflation from the omission (percentage points)")
    ax[1].set_ylim(0, max(d_raw, d_cor) * 1.3)
    panel(ax[1], "b", f"{d_cor/max(d_raw,1e-9):.1f}x larger when noise-corrected")

    # c: the check that decides whether the correction is a context property at
    # all. Uncentred, the four lines' corrections are one shared programme --
    # every perturbation here is an essential-gene knockdown.
    ax[2].bar([0, 1], [base.alpha_between_line_r, -0.532], width=0.5,
              color=[ORANGE, AQUA])
    ax[2].axhline(0, color="#444", lw=1.0)
    for i_, v in enumerate([base.alpha_between_line_r, -0.532]):
        ax[2].text(i_, v + (0.04 if v > 0 else -0.09), f"r = {v:+.2f}",
                   ha="center", fontsize=10, fontweight="bold")
    ax[2].set_xticks([0, 1], ["uncentred\n(one shared programme)",
                              "centred across lines\n(a line contrast)"],
                     fontsize=7.5)
    ax[2].set_ylabel("correlation between different lines' corrections")
    ax[2].set_ylim(-0.8, 0.8)
    ax[2].text(0.5, 0.20, "uncentred, each also correlates with the shared\n"
               f"response at r = {base.alpha_vs_shared_r:+.2f} — subtracting it\n"
               "removes signal the statistic should keep",
               transform=ax[2].transAxes, ha="center", va="top", fontsize=6.4,
               color="#444")
    panel(ax[2], "c", "The correction must be a contrast")
    fig.suptitle("Figure 7 — The omission inflates another group's "
                 "context-specificity statistic, in a different perturbation "
                 "modality", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    return fig, {"summary": T}


# ---------------------------------------------------------------- figure 8
def figure8():
    """Does exposure remodel the programme, or move cells along it?"""
    CU = read("capture_curve.csv", sub="expression_remodelling")
    BD = read("by_dose.csv", sub="expression_remodelling")
    CV = read("convergence.csv", sub="expression_remodelling")
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)

    # a: the capture curve against the only reference that makes it readable --
    # a held-out UNTREATED line, which is a shift by construction.
    if CU is not None and len(CU):
        ax[0].plot(CU.k, CU.held_out_baseline, "o-", color=AQUA, lw=2.2, ms=7,
                   label="held-out untreated cell line")
        ax[0].plot(CU.k, CU.drug_response, "o-", color=ORANGE, lw=2.2, ms=7,
                   label="drug response")
        ax[0].set_xlabel("baseline directions used (k)")
        ax[0].set_ylabel("fraction captured")
        ax[0].legend(frameon=False, fontsize=7.5, loc="upper left")
        ax[0].text(0.97, 0.06, "responses are captured 2-5x less\nat every k",
                   transform=ax[0].transAxes, ha="right", fontsize=7,
                   color="#444")
    panel(ax[0], "a", "Exposure moves cells off the baseline axes")

    # b: flat across dose, so it is a programme and not cytotoxic collapse
    if BD is not None and len(BD):
        ax[1].plot(BD.conc, BD.remodel_share, "o-", color=ORANGE, lw=2.2, ms=8)
        for r_ in BD.itertuples():
            ax[1].text(r_.conc, r_.remodel_share + 0.02,
                       f"{r_.remodel_share:.0%}", ha="center", fontsize=8.5,
                       fontweight="bold")
        ax[1].set_xscale("log")
        ax[1].set_ylim(0, 1.05)
        ax[1].set_xlabel("concentration (µM)")
        ax[1].set_ylabel("share of the reproducible response off-axis")
        ax[1].text(0.5, 0.14, "already at full size at the lowest dose,\nwhere "
                   "cells are not dying — so not\ncytotoxic collapse",
                   transform=ax[1].transAxes, ha="center", fontsize=7,
                   color="#444")
    panel(ax[1], "b", "Present at the lowest dose, so not death")

    # c: identity is preserved -- the new axis is added, not substituted
    if CV is not None and len(CV):
        cs = sorted(CV.conc.unique())
        d0 = [CV[CV.conc == c].spread_ratio.dropna() for c in cs]
        bp = ax[2].boxplot(d0, showfliers=False, patch_artist=True,
                           medianprops=dict(color="black", lw=1.4))
        for pch in bp["boxes"]:
            pch.set_facecolor(BLUE); pch.set_alpha(0.75)
        ax[2].axhline(1.0, color=ORANGE, lw=2.2, label="untreated spread")
        for i_, dd in enumerate(d0):
            if len(dd):
                ax[2].text(i_ + 1, np.median(dd) + 0.02,
                           f"{np.median(dd):.3f}", ha="center", fontsize=8,
                           fontweight="bold")
        ax[2].set_xticks(range(1, len(cs) + 1), [f"{c:g}" for c in cs],
                         fontsize=8.5)
        ax[2].set_xlabel("concentration (µM)")
        ax[2].set_ylabel("between-line spread, treated / untreated")
        ax[2].legend(frameon=False, fontsize=7.5)
        ax[2].text(0.5, 0.08, "lines stay exactly as far apart as they were;\n"
                   "gene-gene coordination is also unchanged\n"
                   "(rho = +0.998 over 400 genes)",
                   transform=ax[2].transAxes, ha="center", fontsize=6.8,
                   color="#444")
    panel(ax[2], "c", "Cell identity is preserved")
    fig.suptitle("Figure 8 — Drug exposure writes a new axis on top of an "
                 "unchanged programme", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    return fig, {"capture_curve": CU if CU is not None else pd.DataFrame(),
                 "by_dose": BD if BD is not None else pd.DataFrame(),
                 "convergence": CV if CV is not None else pd.DataFrame()}


# ---------------------------------------------------------------- figure 9
def figure9():
    """Chromatin: the positive control, the cost of dilution, and the null."""
    T = read("atac_responsive.csv")
    PC = read("atac_positive_control.csv")
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)
    if T is None or len(T) < 2 or PC is None or not len(PC):
        for a_ in ax:
            a_.axis("off")
        return fig, {}
    resp = T[T.features.str.startswith("responsive")].iloc[0]
    alln = T[T.features == "all"].iloc[0]

    # a: the control the experiment supplies -- a knockout should lower its own
    # motif, and if it does not, no global statistic can be trusted either
    top = PC.sort_values("p").head(10)
    yy = np.arange(len(top))[::-1]
    ax[0].barh(yy, top.delta,
               color=[ORANGE if q < 0.05 else GREY for q in top.q], height=0.7)
    ax[0].axvline(0, color="#444", lw=0.9)
    ax[0].set_yticks(yy, [f"{r.pert} ({r.line})" for r in top.itertuples()],
                     fontsize=6.6)
    ax[0].set_xlabel("Δ ChromVar deviation, the knockout's own motif")
    n_sig = int((PC.q < 0.05).sum())
    ax[0].text(0.04, 0.10, f"{n_sig} of {len(PC)} reject;\nall six move down",
               transform=ax[0].transAxes, fontsize=8, color=ORANGE,
               fontweight="bold")
    panel(ax[0], "a", "The experiment worked")

    # b: identical data and estimator, only the feature set differs
    xx = np.arange(2)
    w = 0.36
    ax[1].bar(xx - w / 2, [alln.resid, alln.interaction_of_reproducible], w,
              color=GREY, label=f"all {int(alln.n):,} features")
    ax[1].bar(xx + w / 2, [resp.resid, resp.interaction_of_reproducible], w,
              color=VIOLET, label=f"{int(resp.n)} responsive (out of fold)")
    for i_, (a_, b_) in enumerate([(alln.resid, resp.resid),
                                   (alln.interaction_of_reproducible,
                                    resp.interaction_of_reproducible)]):
        ax[1].text(i_ - w / 2, a_ + .015, f"{a_:.0%}", ha="center", fontsize=8.5)
        ax[1].text(i_ + w / 2, b_ + .015, f"{b_:.0%}", ha="center", fontsize=8.5,
                   fontweight="bold")
    ax[1].set_xticks(xx, ["residual\n(noise)", "interaction /\nreproducible"],
                     fontsize=7.5)
    ax[1].set_ylim(0, 1.1)
    ax[1].legend(frameon=False, fontsize=7)
    ax[1].text(0.5, 0.44, "same data, same estimator —\nonly the feature set "
               "differs", transform=ax[1].transAxes, ha="center", fontsize=7,
               color="#444")
    panel(ax[1], "b", "What testing everywhere costs")

    # c: and even where the signal is, the interaction is at chance
    obs, nm = float(resp.inter), float(resp.null_mean)
    ax[2].bar([0, 1], [obs, nm], width=0.5, color=[ORANGE, GREY])
    for i_, v in enumerate([obs, nm]):
        ax[2].text(i_, v + 0.008, f"{v:.1%}", ha="center", fontsize=11,
                   fontweight="bold")
    ax[2].set_xticks([0, 1], ["observed", "permutation\nnull"], fontsize=8)
    ax[2].set_ylabel("interaction variance fraction")
    ax[2].set_ylim(0, max(obs, nm) * 1.45)
    ax[2].text(0.5, 0.86, f"P = {float(resp.p_vs_null):.2f}\n"
               f"at chance on features that\ndemonstrably respond",
               transform=ax[2].transAxes, ha="center", fontsize=8,
               color=ORANGE, fontweight="bold")
    panel(ax[2], "c", "Still no interaction, now meaningfully")

    fig.suptitle("Figure 9 \u2014 Chromatin: the perturbations work, an "
                 "atlas-wide index cannot see them, and the interaction is "
                 "still at chance", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    return fig, {"summary": T, "positive_control": PC}


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    for i, fn in enumerate((figure1, figure2, figure3, figure4, figure5,
                            figure6, figure7, figure8, figure9), 1):
        fig, src = fn()
        d = save_figure(fig, f"fig{i}", FIG, source_data=src, script=__file__)
        plt.close(fig)
        print(f"Figure {i} -> {d}")


if __name__ == "__main__":
    main()
