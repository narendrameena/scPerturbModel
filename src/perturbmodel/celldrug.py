"""What is the relation between one cell line and one drug?

A single measured response carries three things that are usually collapsed
together, and separating them is the whole problem:

    y(line, drug, dose) = beta(drug, dose)      the drug's average effect
                        + alpha(line)           how sensitive this line is to
                                                everything -- growth rate,
                                                seeding, drug metabolism
                        + gamma(line, drug)     the part that is specific to
                                                THIS pairing
                        + noise

Only ``gamma`` is a cell-drug relation. ``alpha`` is a property of the cell that
has nothing to do with the drug, and ``beta`` is a property of the drug that has
nothing to do with the cell. Every project analysis until now removed ``beta``
and kept ``alpha + gamma`` in the residual, so a line that is simply frail
looked like a line with a specific relation to every compound it was screened
against. Ben-David et al. (Nature 2018) name this directly: their most resistant
MCF7 strains are resistant *in general*, through downregulated drug-metabolism
pathways, not through anything drug-specific.

Two rules make the estimates honest:

  * ``beta`` for a query pair is computed from OTHER lines (leave-one-line-out),
    so the pair cannot contribute to the mean it is measured against.
  * ``alpha`` for a query pair is computed from OTHER compounds
    (leave-one-compound-out), so a genuine relation cannot be absorbed into the
    line's general sensitivity and disappear.

Without the second rule ``alpha`` eats part of ``gamma`` whenever a line has few
compounds, and the interaction is biased toward zero.

``gamma`` is then validated against independent replicate detection plates: the
covariance of gamma between plates X1 and X2 at matched line, compound and dose
is an estimate of the true interaction variance that noise cannot inflate, since
the noise in two plates is independent. A per-pair z uses the same replicate
structure, so a reported relation is one that reproduced, not one that was large
once.

Command line::

    python -m perturbmodel.celldrug --cell ACH-000019 --drug dabrafenib
    python -m perturbmodel.celldrug --cell "MCF7" --top 15
    python -m perturbmodel.celldrug --variance          # global apportionment
"""
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PRISM = ROOT / "data" / "external" / "prism"
MIN_LINES = 30
MIN_CPDS = 20


# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------
def load_prism(prism_dir: Path | None = None):
    """Response cube as (lines, conditions) plus a condition key table.

    Conditions are (compound, dose, replicate detection plate). Keeping the
    replicate axis separate is what later allows the interaction to be validated
    rather than asserted.
    """
    p = Path(prism_dir) if prism_dir else PRISM
    lfc = pd.read_csv(p / "secondary-screen-logfold-change.csv", index_col=0)
    # 8 rows are pool_line_FAILED_STR: unauthenticated cultures that also appear
    # as clean rows. Drop rather than average into the authenticated culture.
    lfc = lfc[[not i.endswith("_FAILED_STR") for i in lfc.index]]
    lfc.index = [re.search(r"(ACH-\d+)", i).group(1) for i in lfc.index]
    lfc = lfc.groupby(level=0).mean()
    ti = pd.read_csv(p / "secondary-screen-replicate-treatment-info.csv",
                     low_memory=False)
    ti = ti[ti.column_name.isin(lfc.columns) & ti.name.notna()].copy()
    ti["rep"] = ti.detection_plate.astype(str).str.extract(r"_(X\d)")[0]
    ti = ti[ti.rep.notna()]
    ti["dose_s"] = ti.dose.round(4).astype(str)
    L = lfc.to_numpy(dtype=np.float32)
    colpos = {c: i for i, c in enumerate(lfc.columns)}
    keys, mats = [], []
    for (cpd, dose, rep), g in ti.groupby(["name", "dose_s", "rep"],
                                          observed=True):
        idx = [colpos[c] for c in g.column_name]
        with np.errstate(invalid="ignore"):
            v = np.nanmean(L[:, idx], axis=1) if len(idx) > 1 else L[:, idx[0]]
        keys.append((cpd, dose, rep))
        mats.append(v)
    K = pd.DataFrame(keys, columns=["compound", "dose", "rep"])
    return np.stack(mats, axis=1), K, np.array(lfc.index), ti


# --------------------------------------------------------------------------
# the decomposition
# --------------------------------------------------------------------------
@dataclass
class Decomposition:
    """Out-of-fold three-way split of a response cube."""
    lines: np.ndarray
    K: pd.DataFrame
    drug_effect: pd.DataFrame        # beta, conditions x 1 (leave-one-line-out)
    line_effect: pd.Series           # alpha, general sensitivity per line
    gamma: dict                      # compound -> Series over lines
    gamma_rep: dict                  # compound -> DataFrame lines x replicate
    alpha_split: pd.DataFrame        # alpha on two disjoint compound halves

    def general_sensitivity(self, line):
        return float(self.line_effect.get(line, np.nan))


def decompose(R, K, lines, min_lines=MIN_LINES, min_cpds=MIN_CPDS, seed=0):
    """Split the cube into drug effect, line general sensitivity, interaction.

    Returns per-replicate interaction so downstream code can test reproducibility
    rather than trusting a single measurement.
    """
    rng = np.random.default_rng(seed)
    nL = R.shape[0]

    # --- beta: leave-one-line-out mean across lines, per (compound, dose, rep)
    with np.errstate(invalid="ignore"):
        colsum = np.nansum(R, axis=0)
        colcnt = np.sum(np.isfinite(R), axis=0)
    ok_col = colcnt > min_lines
    beta = np.full(R.shape, np.nan, dtype=np.float32)
    with np.errstate(invalid="ignore"):
        beta[:, ok_col] = ((colsum[ok_col] - np.nan_to_num(R[:, ok_col]))
                           / (colcnt[ok_col] - 1))
    beta[~np.isfinite(R)] = np.nan
    D = R - beta                       # drug effect removed, alpha + gamma left

    # --- alpha: leave-one-COMPOUND-out mean of D across that line's compounds
    cpd_of_col = K.compound.to_numpy()
    ucpd = pd.unique(cpd_of_col)
    cidx = {c: np.where(cpd_of_col == c)[0] for c in ucpd}
    # per (line, compound) mean of D, then average over compounds excluding self
    per_cpd = np.full((nL, len(ucpd)), np.nan, dtype=np.float32)
    for j, c in enumerate(ucpd):
        with np.errstate(invalid="ignore"):
            per_cpd[:, j] = np.nanmean(D[:, cidx[c]], axis=1)
    with np.errstate(invalid="ignore"):
        rs = np.nansum(per_cpd, axis=1)
        rc = np.sum(np.isfinite(per_cpd), axis=1)
    enough = rc >= min_cpds
    alpha = pd.Series(np.where(enough, rs / np.maximum(rc, 1), np.nan),
                      index=lines, name="general_sensitivity")

    # alpha estimated twice on disjoint compound halves -- the direct test of
    # whether general sensitivity is a reproducible property of the line
    perm = rng.permutation(len(ucpd))
    hA, hB = perm[:len(ucpd) // 2], perm[len(ucpd) // 2:]
    with np.errstate(invalid="ignore"):
        aA = np.nanmean(per_cpd[:, hA], axis=1)
        aB = np.nanmean(per_cpd[:, hB], axis=1)
    alpha_split = pd.DataFrame({"half_A": aA, "half_B": aB}, index=lines)

    # --- alpha PER REPLICATE PLATE. The log-fold-change is taken against
    # control wells on the same detection plate, so every condition a line has
    # on plate X1 carries that plate's control noise. An alpha averaged over
    # plates carries a share of each, and the cross-terms then subtract part of
    # that nuisance from the replicate covariance below rather than leaving it
    # to cancel. Estimated within a plate it removes exactly the term its own
    # residual carries, and the two members of a cross-plate pair are corrected
    # from independent plates.
    reps = K.rep.to_numpy()
    alpha_by_rep = {}
    for rp in pd.unique(reps):
        rc_cols = np.where(reps == rp)[0]
        cpd_r = cpd_of_col[rc_cols]
        pcr = np.full((nL, len(ucpd)), np.nan, dtype=np.float32)
        for j, c in enumerate(ucpd):
            sel = rc_cols[cpd_r == c]
            if len(sel):
                with np.errstate(invalid="ignore"):
                    pcr[:, j] = np.nanmean(D[:, sel], axis=1)
        with np.errstate(invalid="ignore"):
            rs_r = np.nansum(pcr, axis=1)
            rc_r = np.sum(np.isfinite(pcr), axis=1)
        a_r = np.full((nL, len(ucpd)), np.nan, dtype=np.float32)
        good = rc_r >= 5
        with np.errstate(invalid="ignore"):
            a_r[good] = ((rs_r[good, None] - np.nan_to_num(pcr[good]))
                         / np.maximum(rc_r[good, None] - 1, 1))
        alpha_by_rep[rp] = a_r

    # --- gamma: what is left once BOTH main effects are out
    gamma, gamma_rep = {}, {}
    for j, c in enumerate(ucpd):
        cols = cidx[c]
        sub = K.iloc[cols]
        per_rep = {}
        for rep, gr in sub.groupby("rep", observed=True):
            pos = [list(cols).index(i) for i in gr.index]
            a_r = alpha_by_rep[rep][:, j]
            g = D[:, cols][:, pos] - a_r[:, None]
            with np.errstate(invalid="ignore"):
                per_rep[rep] = np.nanmean(g, axis=1)
        Gr = pd.DataFrame(per_rep, index=lines).dropna(how="all")
        if Gr.shape[1] >= 2 and Gr.notna().all(axis=1).sum() >= min_lines:
            gamma_rep[c] = Gr
            gamma[c] = Gr.mean(axis=1).dropna()
    return Decomposition(lines, K, pd.DataFrame({"beta_mean": np.nanmean(beta, 0)}),
                         alpha, gamma, gamma_rep, alpha_split)


# --------------------------------------------------------------------------
# the canonical corrected residual, shared by every downstream analysis
# --------------------------------------------------------------------------
CACHE = ROOT / "data" / "processed" / "prism_gamma.pkl"


def prism_gamma(cache: Path | None = None, rebuild: bool = False):
    """Interaction residuals per compound, with general sensitivity removed.

    This is the drop-in replacement for the old project convention, in which the
    residual was the response minus the leave-one-line-out compound mean and so
    still contained the line's general sensitivity. Every analysis that asks a
    question about a cell-drug relation should take its residual from here, so
    that the correction lives in one place and cannot drift between scripts.

    Returns ``(gamma, ti, dec)`` where ``gamma`` maps compound -> Series over
    cell lines, matching the shape the old ``residuals()`` returned.

    The decomposition costs several minutes over the full PRISM cube, so the
    result is cached; pass ``rebuild=True`` after changing the estimator.
    """
    import pickle

    c = Path(cache) if cache else CACHE
    if c.exists() and not rebuild:
        with open(c, "rb") as fh:
            d = pickle.load(fh)
        return d["gamma"], d["ti"], d["dec"]
    R, K, lines, ti = load_prism()
    dec = decompose(R, K, lines)
    c.parent.mkdir(parents=True, exist_ok=True)
    with open(c, "wb") as fh:
        pickle.dump({"gamma": dec.gamma, "ti": ti, "dec": dec}, fh)
    return dec.gamma, ti, dec


def general_sensitivity_by_rep(R, K, rows=None, min_cpds=5, rep_col="rep"):
    """Per-replicate-plate general sensitivity, leave-one-compound-out.

    The correction a residual gets must be built from the SAME plate that
    residual came from. A log-fold-change is taken against control wells on its
    own detection plate, so every condition a line has on that plate carries the
    plate's control noise; a correction pooled over plates carries a share of
    each and subtracts it from the wrong side of a cross-plate replicate pair.
    In simulation that removed var(control)/2 from the interaction covariance
    and drove the estimate to exactly zero.

    Returns ``(by_rep, order)`` where ``by_rep[r][:, j]`` is the correction for
    plate ``r`` and compound ``order[j]``. Use ``general_sensitivity_loo`` only
    where the residual is already pooled across plates.
    """
    order = [c for c, _ in K.groupby("compound", observed=True)]
    reps = K[rep_col].to_numpy()
    out = {}
    for rp in pd.unique(reps):
        sub = K[reps == rp]
        A, sub_order = general_sensitivity_loo(R, sub, rows=rows,
                                               min_cpds=min_cpds)
        pos = {c: j for j, c in enumerate(sub_order)}
        n = A.shape[0]
        full = np.full((n, len(order)), np.nan, dtype=np.float32)
        for j, c in enumerate(order):
            if c in pos:
                full[:, j] = A[:, pos[c]]
        out[rp] = full
    return out, order


def general_sensitivity_loo(R, K, rows=None, min_cpds=5):
    """Leave-one-compound-out general sensitivity, for cube-based analyses.

    Returns ``(alpha, order)`` where ``alpha[:, j]`` is each line's mean response
    across every compound EXCEPT ``order[j]``, on the compound-centred scale. Any
    analysis that residualises a response cube against a compound mean should
    subtract the matching column before calling what remains an interaction.

    Leaving the compound out matters: with it in, a line screened against few
    compounds has its own interaction absorbed into its general sensitivity and
    the interaction is biased toward zero.
    """
    idx = np.arange(R.shape[0]) if rows is None else np.asarray(rows)
    groups = list(K.groupby("compound", observed=True))
    with np.errstate(invalid="ignore"):
        per_cpd = np.stack([np.nanmean(R[np.ix_(idx, g.index.to_numpy())],
                                       axis=1) for _, g in groups], axis=1)
        # centre each compound so alpha is a deviation from the compound mean
        # and not a mixture of that with the compounds' own average potency;
        # otherwise alpha carries a term identical for every line, which shifts
        # the residual rather than describing the line
        per_cpd = per_cpd - np.nanmean(per_cpd, axis=0, keepdims=True)
        csum = np.nansum(per_cpd, axis=1)
        ccnt = np.sum(np.isfinite(per_cpd), axis=1)
    with np.errstate(invalid="ignore"):
        alpha = ((csum[:, None] - np.nan_to_num(per_cpd))
                 / np.maximum(ccnt[:, None] - np.isfinite(per_cpd), 1))
    alpha[ccnt < min_cpds] = np.nan
    return alpha, [c for c, _ in groups]


def remove_line_effect(res_by_cpd, min_cpds=MIN_CPDS):
    """Strip general sensitivity from residuals that already exist.

    For analyses that build their residuals from a source other than the PRISM
    cube -- Tahoe, LINCS, GDSC -- and so cannot call ``prism_gamma``. Takes the
    old-style ``{compound: Series over contexts}`` and returns the same
    structure with each context's across-compound mean removed, computed
    leave-one-compound-out so a genuine interaction is not absorbed.
    """
    M = pd.DataFrame(res_by_cpd)
    cnt = M.notna().sum(axis=1)
    tot = M.sum(axis=1, skipna=True)
    keep = cnt >= min_cpds
    out = {}
    for c in M.columns:
        col = M[c]
        loo = (tot - col.fillna(0.0)) / np.maximum(cnt - col.notna(), 1)
        g = (col - loo)[keep & col.notna()]
        if len(g):
            out[c] = g
    return out


# --------------------------------------------------------------------------
# variance apportionment, replicate-validated
# --------------------------------------------------------------------------
def apportion(dec: Decomposition):
    """How much of the response is drug, cell property, and cell-drug relation.

    Each component is measured as a covariance between two independent estimates
    so that noise, which is uncorrelated between them, contributes zero:

      var(alpha)  cov of the line effect between disjoint compound halves
      var(gamma)  cov of the interaction between independent replicate plates
      var(beta)   variance of the drug effect across conditions

    A variance estimated as a covariance can come out negative when the true
    value is near zero; that is reported as-is rather than clipped, because
    clipping is what turns a null result into a positive one.
    """
    s = dec.alpha_split.dropna()
    v_alpha = float(np.cov(s.half_A, s.half_B)[0, 1])
    covs, weights = [], []
    for c, Gr in dec.gamma_rep.items():
        cols = list(Gr.columns)
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                d = Gr[[cols[i], cols[j]]].dropna()
                if len(d) >= MIN_LINES:
                    covs.append(float(np.mean(d.iloc[:, 0] * d.iloc[:, 1])))
                    weights.append(len(d))
    v_gamma = float(np.average(covs, weights=weights)) if covs else np.nan
    v_beta = float(np.nanvar(dec.drug_effect.beta_mean))
    tot = v_beta + max(v_alpha, 0) + max(v_gamma, 0)
    return {"var_drug": v_beta, "var_cell_property": v_alpha,
            "var_cell_drug_relation": v_gamma,
            "share_drug": v_beta / tot, "share_cell_property": max(v_alpha, 0) / tot,
            "share_relation": max(v_gamma, 0) / tot,
            "n_compounds": len(dec.gamma_rep), "n_lines": int(len(s)),
            "gamma_covs": covs}


# --------------------------------------------------------------------------
# the per-pair query
# --------------------------------------------------------------------------
def relation(dec: Decomposition, line: str, drug: str):
    """The relation between one cell line and one drug, with a reproducibility z.

    ``z`` compares the interaction against the spread of the same quantity
    measured on independent replicate plates for that pair, so it answers "did
    this pairing behave the same way twice", not "was it far from the mean once".
    """
    if drug not in dec.gamma_rep:
        return {"error": f"no validated data for compound {drug!r}"}
    Gr = dec.gamma_rep[drug]
    if line not in Gr.index:
        return {"error": f"line {line!r} not screened against {drug!r}"}
    vals = Gr.loc[line].dropna()
    if len(vals) < 2:
        return {"error": f"only {len(vals)} replicate(s) for this pair"}
    g = float(vals.mean())
    se = float(vals.std(ddof=1) / np.sqrt(len(vals)))
    pool = dec.gamma[drug]
    z_vs_panel = (g - float(pool.mean())) / max(float(pool.std()), 1e-9)
    return {"line": line, "drug": drug,
            "interaction": g,
            "replicate_se": se,
            "n_replicates": int(len(vals)),
            "z_reproducible": g / se if se > 0 else np.nan,
            "z_vs_panel": z_vs_panel,
            "general_sensitivity": dec.general_sensitivity(line),
            "direction": "more sensitive than the drug's average, beyond this "
                         "line's general sensitivity" if g < 0 else
                         "more resistant than the drug's average, beyond this "
                         "line's general sensitivity",
            "panel_sd": float(pool.std())}


def _resolve(name, lines, ci):
    if name in set(lines):
        return name
    m = ci[ci.ccle_name.astype(str).str.upper().str.startswith(
        str(name).upper() + "_")]
    if len(m):
        return m.depmap_id.iloc[0]
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cell")
    ap.add_argument("--drug")
    ap.add_argument("--top", type=int, default=0,
                    help="list the N strongest relations for --cell or --drug")
    ap.add_argument("--variance", action="store_true",
                    help="global apportionment into drug / cell / relation")
    ap.add_argument("--prism-dir")
    a = ap.parse_args(argv)

    R, K, lines, ti = load_prism(Path(a.prism_dir) if a.prism_dir else None)
    print(f"loaded {R.shape[0]} lines x {R.shape[1]} conditions "
          f"({K.compound.nunique()} compounds)", flush=True)
    dec = decompose(R, K, lines)
    print(f"decomposed: {len(dec.gamma_rep)} compounds with >=2 replicate "
          f"plates", flush=True)

    if a.variance or not (a.cell or a.drug):
        v = apportion(dec)
        s = dec.alpha_split.dropna()
        r = float(np.corrcoef(s.half_A, s.half_B)[0, 1])
        print("\nWhat the response is made of")
        print(f"  drug effect (same in every line)      "
              f"{v['share_drug']:6.1%}   var {v['var_drug']:.4f}")
        print(f"  cell property (general sensitivity)   "
              f"{v['share_cell_property']:6.1%}   var "
              f"{v['var_cell_property']:.4f}")
        print(f"  cell-drug relation (interaction)      "
              f"{v['share_relation']:6.1%}   var "
              f"{v['var_cell_drug_relation']:.4f}")
        print(f"\n  general sensitivity reproduces across disjoint compound "
              f"halves at r = {r:.3f} ({v['n_lines']} lines)")
        print("  -- it is a real property of the line, and it must be removed "
              "before\n     anything is called a cell-drug relation.")

    ci = pd.read_csv(PRISM / "secondary-screen-cell-line-info.csv")
    cell = _resolve(a.cell, lines, ci) if a.cell else None
    if a.cell and cell is None:
        print(f"\ncould not resolve cell {a.cell!r}")
        return
    if cell and a.drug:
        r = relation(dec, cell, a.drug)
        print(f"\n{a.cell} ({cell})  x  {a.drug}")
        for k, val in r.items():
            print(f"  {k:22s} {val:.4f}" if isinstance(val, float)
                  else f"  {k:22s} {val}")
    elif cell and a.top:
        rows = []
        for c in dec.gamma_rep:
            rr = relation(dec, cell, c)
            if "error" not in rr:
                rows.append(rr)
        T = pd.DataFrame(rows)
        T["abs_z"] = T.z_reproducible.abs()
        T = T.sort_values("abs_z", ascending=False).head(a.top)
        print(f"\nstrongest reproducible relations for {a.cell} ({cell}); "
              f"general sensitivity {dec.general_sensitivity(cell):+.3f}")
        print(T[["drug", "interaction", "z_reproducible", "z_vs_panel",
                 "n_replicates"]].round(3).to_string(index=False))
    elif a.drug and a.top:
        Gr = dec.gamma_rep.get(a.drug)
        if Gr is None:
            print(f"no validated data for {a.drug!r}")
            return
        rows = [relation(dec, ln, a.drug) for ln in Gr.index]
        T = pd.DataFrame([r for r in rows if "error" not in r])
        T["abs_z"] = T.z_reproducible.abs()
        T = T.sort_values("abs_z", ascending=False).head(a.top)
        print(f"\nlines with the strongest reproducible relation to {a.drug}")
        print(T[["line", "interaction", "z_reproducible", "general_sensitivity",
                 "n_replicates"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()


def remove_line_effect_profiles(out, min_cpds=3):
    """Strip each line's general response from (line, compound) -> profile.

    The profile analogue of ``remove_line_effect``, for transcriptional arms
    where the residual is a gene vector rather than a scalar. For each line the
    mean of its profiles across the OTHER compounds is subtracted, so that what
    remains is specific to the pairing rather than to how that culture responds
    to being perturbed at all.

    Lines with fewer than ``min_cpds`` compounds are dropped: their general
    response cannot be estimated without absorbing the very interaction being
    measured, and keeping them with a zero correction would silently mix two
    different quantities in one table.
    """
    by_line = {}
    for (ln, cpd), v in out.items():
        by_line.setdefault(ln, {})[cpd] = v
    res = {}
    for ln, d in by_line.items():
        if len(d) < min_cpds:
            continue
        cpds = list(d)
        V = np.stack([d[c] for c in cpds])
        tot, m = V.sum(0), len(V)
        for j, c in enumerate(cpds):
            res[(ln, c)] = V[j] - (tot - V[j]) / (m - 1)
    return res
