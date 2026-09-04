"""A context that responds to everything must not be reported as interaction.

This is the defect these tests exist to prevent, found in August 2026: every
residual in the project removed the compound main effect but not each cell
line's general sensitivity. That term reproduces across disjoint compound halves
at r = 0.989 in PRISM and is shared between replicate detection plates exactly
as a genuine interaction is, so it landed in the pair covariance and was counted
as context-dependence, inflating it by roughly half.

The tests below plant a known general-sensitivity term and check that the
estimators do not pay for it.
"""
import numpy as np
import pandas as pd
import pytest

from perturbmodel.celldrug import (apportion, decompose,
                                   general_sensitivity_loo,
                                   remove_line_effect,
                                   remove_line_effect_profiles)


def _cube(n_line=60, n_cpd=40, n_dose=2, n_rep=2, s_ctx=0.0, s_int=0.5,
          s_noise=0.3, s_plate=0.0, seed=0):
    """Response cube with known drug, line and interaction components.

    ``s_plate`` is a per-(line, plate) offset, the analogue of control-well noise
    in a log-fold-change: every condition a line has on one detection plate
    carries it. It is what makes the correction's own grouping matter -- a
    general-sensitivity estimate pooled over plates carries a share of each
    plate's offset and subtracts it from the wrong side of the replicate pair.
    """
    rng = np.random.default_rng(seed)
    beta = rng.normal(0, 1.0, (n_cpd, n_dose))
    alpha = rng.normal(0, s_ctx, n_line)
    gamma = rng.normal(0, s_int, (n_line, n_cpd))
    plate = rng.normal(0, s_plate, (n_line, n_rep))
    keys, cols = [], []
    for c in range(n_cpd):
        for d in range(n_dose):
            for r in range(n_rep):
                keys.append((f"cpd{c}", str(d), f"X{r}"))
                cols.append(beta[c, d] + alpha + gamma[:, c] + plate[:, r]
                            + rng.normal(0, s_noise, n_line))
    K = pd.DataFrame(keys, columns=["compound", "dose", "rep"])
    R = np.stack(cols, axis=1).astype(np.float32)
    lines = np.array([f"ACH-{i:06d}" for i in range(n_line)])
    return R, K, lines, alpha, gamma


def test_general_sensitivity_is_recovered():
    """alpha is estimated, and it tracks the planted values."""
    R, K, lines, alpha, _ = _cube(s_ctx=1.0, seed=1)
    dec = decompose(R, K, lines, min_lines=10, min_cpds=5)
    est = dec.line_effect.reindex(lines).to_numpy()
    ok = np.isfinite(est)
    assert ok.sum() > 40
    assert np.corrcoef(est[ok], alpha[ok])[0, 1] > 0.9


def test_interaction_does_not_absorb_general_sensitivity():
    """The reported interaction must not grow when only alpha grows.

    This is the exact failure: with the line term left in the residual, its
    variance is added to the pair covariance, so a panel of frail lines reads as
    a panel of context-specific ones.
    """
    quiet = apportion(decompose(*_cube(s_ctx=0.0, seed=2)[:3],
                                min_lines=10, min_cpds=5))
    loud = apportion(decompose(*_cube(s_ctx=1.5, seed=2)[:3],
                               min_lines=10, min_cpds=5))
    # a 1.5-SD general-sensitivity term is planted in `loud` and none in
    # `quiet`; the interaction variance must be essentially unchanged
    assert loud["var_cell_drug_relation"] == pytest.approx(
        quiet["var_cell_drug_relation"], rel=0.25)
    # and it must be picked up as a cell property instead
    assert loud["var_cell_property"] > 5 * max(quiet["var_cell_property"], 1e-4)


def test_interaction_variance_is_close_to_truth():
    R, K, lines, _, gamma = _cube(s_ctx=1.0, s_int=0.5, seed=3)
    v = apportion(decompose(R, K, lines, min_lines=10, min_cpds=5))
    assert v["var_cell_drug_relation"] == pytest.approx(gamma.var(), rel=0.3)


def test_loo_helper_matches_decompose():
    """general_sensitivity_loo is what cube-based scripts subtract."""
    R, K, lines, alpha, _ = _cube(s_ctx=1.0, seed=4)
    A, order = general_sensitivity_loo(R, K)
    assert A.shape == (R.shape[0], K.compound.nunique())
    assert len(order) == K.compound.nunique()
    m = np.nanmean(A, axis=1)
    assert np.corrcoef(m, alpha)[0, 1] > 0.9


def test_leave_one_compound_out_excludes_self():
    """A compound must not contribute to its own correction."""
    R, K, lines, _, _ = _cube(n_cpd=10, s_ctx=1.0, seed=5)
    A, order = general_sensitivity_loo(R, K)
    j = order.index("cpd3")
    R2 = R.copy()
    cols = K.index[K.compound == "cpd3"].to_numpy()
    R2[:, cols] += 100.0                       # enormous shift, one compound
    A2, _ = general_sensitivity_loo(R2, K)
    assert np.allclose(A[:, j], A2[:, j], atol=1e-3)


def test_remove_line_effect_scalar_and_profile():
    rng = np.random.default_rng(6)
    lines = [f"L{i}" for i in range(30)]
    a = dict(zip(lines, rng.normal(0, 1, 30)))
    scal = {f"c{k}": pd.Series({l: a[l] + rng.normal(0, .1) for l in lines})
            for k in range(25)}
    out = remove_line_effect(scal)
    assert max(abs(float(np.mean(v))) for v in out.values()) < 0.5
    prof = {(l, f"c{k}"): a[l] * np.ones(5) + rng.normal(0, .1, 5)
            for l in lines for k in range(6)}
    outp = remove_line_effect_profiles(prof)
    assert len(outp) == len(prof)
    assert np.abs(np.stack(list(outp.values()))).mean() < 0.4


@pytest.mark.parametrize("s_plate", [0.0, 0.8, 1.5])
def test_per_plate_control_noise_does_not_destroy_the_interaction(s_plate):
    """A per-(line, plate) offset must not eat the interaction.

    The regression this guards: the general-sensitivity correction was first
    estimated by pooling a line's compounds across BOTH plates, so it carried
    half of each plate's control noise. The cross-terms of the replicate
    covariance then subtracted var(control)/2, which drove the estimated
    interaction to exactly zero in simulation while looking perfectly reasonable
    in code review.
    """
    R, K, lines, _, gamma = _cube(s_ctx=1.0, s_int=0.5, s_plate=s_plate, seed=7)
    v = apportion(decompose(R, K, lines, min_lines=10, min_cpds=5))
    assert v["var_cell_drug_relation"] > 0.5 * gamma.var()
    assert v["var_cell_drug_relation"] == pytest.approx(gamma.var(), rel=0.35)


def test_anova_components_recovers_known_variance():
    """The published estimator must recover components it was given.

    Henderson/ANOVA method of moments on a balanced two-way crossed random
    design (Searle, Casella & McCulloch 1992) -- the balanced-design limit of
    variancePartition, used to check this project's bespoke covariance
    estimator against an independent published route. A context term estimated
    from only 3 levels is poorly determined and is not asserted here; the
    interaction, which is the quantity in question, is.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    from published_methods_check import anova_components
    rng = np.random.default_rng(0)
    a, b, n, F = 3, 41, 2, 400
    va, vb, vab, ve = 1.0, 4.0, 0.25, 0.09
    Y = (rng.normal(0, np.sqrt(va), (a, 1, 1, F))
         + rng.normal(0, np.sqrt(vb), (1, b, 1, F))
         + rng.normal(0, np.sqrt(vab), (a, b, 1, F))
         + rng.normal(0, np.sqrt(ve), (a, b, n, F)))
    VP, neg = anova_components(Y, a, b, n)
    tot = va + vb + vab + ve
    assert VP["inter"].median() == pytest.approx(vab / tot, abs=0.02)
    assert VP["resid"].median() == pytest.approx(ve / tot, abs=0.02)
    r = VP["inter"] / (VP["ctx"] + VP["pert"] + VP["inter"])
    assert float(r.median()) == pytest.approx(vab / (va + vb + vab), abs=0.02)


def test_edistance_is_zero_for_identical_distributions():
    """E-distance (Peidli et al. 2024) must vanish when the two samples match.

    The within-group term E||x - x'|| is over DISTINCT pairs; averaging the full
    n x n distance matrix includes the zero diagonal and biases that term low,
    inflating the statistic. In 2,174 dimensions at n = 400 that bias is ~0.33 --
    larger than most real effects here, so it would have made every perturbation
    look significant.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    from published_methods_check import edistance, etest
    rng = np.random.default_rng(0)
    X = rng.normal(0, 1, (300, 500))
    Z = rng.normal(0, 1, (300, 500))
    assert abs(edistance(X, Z)) < 0.05
    # and it must grow with a real shift, and its test must be calibrated
    Y = rng.normal(0.15, 1, (300, 500))
    assert edistance(X, Y) > 5 * abs(edistance(X, Z))
    _, p_null = etest(X, Z, n_perm=60)
    assert p_null > 0.05
