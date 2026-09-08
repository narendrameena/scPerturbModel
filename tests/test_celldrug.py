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


def test_trace_estimator_is_unbiased_under_concentration():
    """The replicate-covariance trace must not depend on how concentrated the
    interaction is.

    This guards a claim that was ASSERTED and then refuted by simulation: that
    the published-style index collapses when the interaction occupies few
    directions. It does not, and a proposed spectral replacement was rejected on
    this evidence (RESULTS.md sec.38). If a future change makes the trace
    concentration-dependent, that is a regression.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    from spectrum_benchmark import simulate, trace_index
    for rank in (2000, 100, 4):
        A, B, truth = simulate(rank=rank, seed=0)
        assert trace_index(A, B) == pytest.approx(truth, rel=0.05)


def test_naive_variance_index_is_badly_wrong_at_the_null():
    """Residual variance reports a large interaction from data containing none.

    Kept as an executable statement of why replicate covariance is required:
    the variance-based alternative returns ~0.5 where the truth is 0.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    from spectrum_benchmark import variance_index, trace_index
    rng = np.random.default_rng(0)
    A = rng.normal(0, 1.0, (300, 500))
    B = rng.normal(0, 1.0, (300, 500))
    assert variance_index(A, B) > 0.4
    assert abs(trace_index(A, B)) < 0.02


def test_sparse_detection_beats_pooled_when_concentrated():
    """Higher Criticism must find sparse interaction the pooled test misses.

    The claim of RESULTS.md sec.39, as an executable check. It is a statement
    about POWER RATES, so it is checked across seeds: the pooled test retains
    about 12% power in this regime rather than 0%, and a first version of this
    test asserted that it must fail on one particular seed -- which is asserting
    the outcome of a coin flip and duly failed on seed 1.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    from sparse_benchmark import simulate_sparse
    from perturbmodel.sparse_interaction import detect
    hc, pooled, localised = 0, 0, 0
    seeds = range(4)
    for seed in seeds:
        A, B, truth = simulate_sparse(k=4, seed=seed)
        r = detect(A, B, n_perm=200, seed=seed)
        hc += r.hc_p < 0.05
        pooled += r.trace_p < 0.05
        hits = set(np.where((r.table.q < 0.05).to_numpy())[0])
        localised += len(hits & truth) >= 1
    assert hc == len(list(seeds))        # sparse test fires every time
    assert pooled < hc                   # pooled test does not keep up
    assert localised == len(list(seeds))  # and every run localises a true pair


def test_sparse_detection_is_clean_at_the_null():
    """Neither selection procedure may fire when there is nothing there."""
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    from sparse_benchmark import simulate_sparse
    from perturbmodel.sparse_interaction import detect
    A, B, _ = simulate_sparse(k=0, total=0.0, seed=7)
    r = detect(A, B, n_perm=200, seed=7)
    assert r.hc_p > 0.05
    assert r.n_sig == 0


def test_permutation_null_crosses_blocks():
    """The null pairing must not systematically re-pair within a block.

    The bug this guards: conditions arrive sorted by cell line, so a cyclic
    shift paired 99.3% of null pairs within the same line against 2.1% at
    chance, putting leftover context structure inside the null (RESULTS.md
    sec.40).
    """
    from perturbmodel.sparse_interaction import _derangement
    rng = np.random.default_rng(0)
    n = 600
    block = np.repeat(np.arange(20), n // 20)      # sorted, as the data are
    idx = _derangement(n, rng, block)
    assert (idx == np.arange(n)).sum() == 0        # no fixed points
    assert (block == block[idx]).mean() < 0.02     # and it crosses blocks


def test_design_calculator_ranks_atlases_correctly():
    """The design calculation must reproduce which atlases resolved anything.

    RESULTS.md sec.46: fed only each atlas's context, perturbation and replicate
    counts, the calculator must call Tahoe and Spear-ATAC unresolvable and LINCS,
    OP3 and sci-Plex resolvable -- the outcomes this project spent weeks
    establishing empirically.
    """
    from perturbmodel.design import audit_design, n_pairs
    # ONE shared noise constant for all five -- a per-atlas value would let
    # five free parameters fit five binary outcomes.
    cases = [("Tahoe-100M", 48, 95, 2, 0.005, False),
             ("LINCS phase 1", 71, 831, 3, 0.57, True),
             ("OP3", 6, 147, 3, 0.331, True),
             ("sci-Plex 3", 3, 189, 2, 0.302, True),
             ("Spear-ATAC", 3, 41, 5, 0.014, False)]
    for snr in (0.10, 0.20, 0.30):        # 5/5 must hold across the range
        for name, c, p_, r, obs, truth in cases:
            a = audit_design(name, c, p_, r, snr, observed_share=obs)
            assert a["resolvable"] is truth, f"{name} at snr={snr}"
    # cells are absent by construction: only replicate pairs enter
    assert n_pairs(10, 10, 1) == 0
    assert n_pairs(10, 10, 3) == 3 * n_pairs(10, 10, 2)


def test_context_embedding_rule():
    """d >= 0.05 * n_contexts should capture essentially all the interaction."""
    from perturbmodel.design import captured_fraction
    for n_ctx in (20, 100, 200, 1000):
        d = int(np.ceil(0.05 * n_ctx)) + 1
        assert captured_fraction(d, n_ctx) > 0.95
    # and a d chosen well below the line does not
    assert captured_fraction(5, 1000) < 0.2


def test_replicated_fraction_changes_the_verdict():
    """Quoting total perturbations instead of replicated ones flips the call.

    RESULTS.md sec.48. Tahoe profiles ~1100 compounds but replicates 13.5% of
    conditions, and only the replicated ones contribute cross-replicate pairs.
    The CLI must not let the total be mistaken for the effective count: at 1100
    the design looks comfortable, at the true 148 it is underpowered against
    Tahoe's own measured share of 0.005.
    """
    from perturbmodel.design import min_detectable_share
    naive = min_detectable_share(48, 1100, 2, 0.20)
    honest = min_detectable_share(48, int(round(1100 * 0.135)), 2, 0.20)
    assert naive < 0.005 < honest, (naive, honest)
    # the inflation is the square root of the perturbation ratio
    assert 2.5 < honest / naive < 2.9


def test_one_replicate_is_never_optimal_under_any_budget():
    """sec.48's budget advice must not depend on the noise-curve parameters.

    The two parameters of the saturating snr curve are assumptions. The
    qualitative claim -- that spending a fixed cell budget on a single replicate
    is never optimal -- has to survive sweeping them, because one replicate
    yields zero pairs however the curve is modelled.
    """
    from perturbmodel.design import min_detectable_share

    def best_reps(cells, n_ctx, n_pert, snr_max, snr_half, max_r=12):
        out = None
        for r in range(1, max_r + 1):
            per = int(cells / (n_ctx * n_pert * r))
            if per < 50:
                break
            snr = snr_max / (1.0 + snr_half / max(per, 1))
            m = min_detectable_share(n_ctx, n_pert, r, snr)
            if out is None or m < out[1]:
                out = (r, m)
        return out[0]

    for snr_half in (50, 100, 300, 1000, 2000):
        for snr_max in (0.1, 0.35, 0.8):
            r = best_reps(95_600_000, 48, 1100, snr_max, snr_half)
            assert r >= 2, (snr_half, snr_max, r)


def test_missing_levels_are_not_counted_as_contexts():
    """sec.48's reader must not count '' or 'None' as a real level.

    Unlabelled cells become '' when a categorical code is -1 and some datasets
    write the literal string 'None'. Counting either inflates the design:
    sci-Plex 3 reads as 4 cell lines and 3 replicates when it has 3 and 2, which
    would overstate its pair count by half.
    """
    from perturbmodel import atlas_meta as m
    col = pd.Series(["MCF7", "A549", "K562", "", "None", "nan", "MCF7"])
    assert m.nlev(col).nunique() == 3
    assert m.pick(["cell_line", "perturbation"], m.CTX) == "cell_line"
    # nperts counts perturbations per cell and must never be read as a replicate
    assert "nperts" not in m.REP


def test_frozen_prereg_artefacts_are_unmodified():
    """The pre-registration is worthless if its frozen files can drift.

    docs/PREREGISTRATION.md records a SHA-256 for each artefact the registered
    analysis depends on. This recomputes them. A failure means either the
    registration must be re-stated or the edit reverted -- it must not pass
    silently.
    """
    import hashlib
    import re
    from pathlib import Path
    root = Path(__file__).resolve().parent.parent
    doc = (root / "docs" / "PREREGISTRATION.md").read_text()
    rows = re.findall(r"\| `([^`]+\.py)` \| `([0-9a-f]{64})` \|", doc)
    assert len(rows) >= 3, "frozen-artefact table missing or malformed"
    for rel, want in rows:
        got = hashlib.sha256((root / rel).read_bytes()).hexdigest()
        assert got == want, (
            f"{rel} changed since pre-registration.\n"
            f"  registered {want}\n  now        {got}")


def test_gem_groups_are_not_independent_replicates():
    """A capture split must not be counted as a replicate.

    RESULTS.md sec.48: Replogle's `batch` is the 10x gem group -- one transduced
    pool distributed across captures -- so two cells with the same guide in
    different gem groups share transduction, culture and selection and differ
    only in the emulsion. Counting them inflates the pair count and would have
    promoted a dataset that cannot support the estimate. The metadata cannot
    distinguish the two cases, so the reported replicate count is an UPPER bound
    and every downstream count must be stated as such.
    """
    from perturbmodel.design import min_detectable_share, n_pairs
    # 2 contexts x 2056 shared perturbations, scored both ways
    as_if_replicated = n_pairs(2, 2056, 40)
    honest = n_pairs(2, 2056, 1)
    assert honest == 0, "one usable replicate must yield no pairs"
    assert as_if_replicated > 3_000_000
    # and the floor difference is the whole verdict
    assert min_detectable_share(2, 2056, 40, 0.20) < 0.001
    assert min_detectable_share(2, 2056, 1, 0.20) == 1.0


def test_five_atlas_table_matches_the_code():
    """The published five-atlas table must not drift from the calculation.

    It did: the manuscript carried floors computed under the earlier per-atlas
    noise constants (Tahoe 0.0079, sci-Plex 0.0109) long after the estimator moved
    to a single shared snr = 0.20, which changes every floor and changes Tahoe's
    prescription from three replicates to six. The verdicts happened to survive,
    which is exactly why nothing caught it. This pins the numbers that appear in
    MANUSCRIPT_DESIGN.md so a future change to the constant fails loudly here.
    """
    from perturbmodel.design import min_detectable_share, n_pairs
    SNR = 0.20
    published = [
        # atlas,          ctx, pert, rep,   pairs,   floor, observed, resolvable
        ("Tahoe-100M",     48,   95,   2,   4_560, 0.0169,    0.005, False),
        ("LINCS phase 1",  71,  831,   3, 177_003, 0.0027,    0.570, True),
        ("OP3",             6,  147,   3,   2_646, 0.0222,    0.331, True),
        ("sci-Plex 3",      3,  189,   2,     567, 0.0479,    0.302, True),
        ("Spear-ATAC",      3,   41,   5,   1_230, 0.0325,    0.014, False),
    ]
    for name, c, p_, r, pairs, floor, obs, resolvable in published:
        assert n_pairs(c, p_, r) == pairs, name
        got = min_detectable_share(c, p_, r, SNR)
        assert abs(got - floor) < 5e-5, f"{name}: table {floor}, code {got:.4f}"
        assert (obs > got) is resolvable, name


def test_no_inner_loop_rebinds_an_enclosing_loop_variable():
    """Guard the bug class that silently broke the cross-laboratory arm.

    In `cross_lab_transcription.load_lincs` an inner loop over cell lines was
    written `for k, ...` while the enclosing loop over compounds was also `k`.
    The compound key was therefore overwritten by a cell line, and every profile
    was stored under (line, last-line-seen) instead of (line, compound). The two
    datasets then shared no keys, so every cross-laboratory comparison returned
    zero pairs -- silently, with exit code 0 and a plausible-looking figure.

    The pattern is only dangerous when the shadowed name is read again after the
    inner loop, which is what this checks. Validated against the pre-fix source:
    it flags that file and no other in the tree.
    """
    import ast
    from pathlib import Path

    def targets(node):
        return {n.id for n in ast.walk(node.target) if isinstance(n, ast.Name)}

    def read_after(outer, name, inner):
        seen = False
        for st in outer.body:
            if st is inner:
                seen = True
                continue
            if not seen:
                continue
            for n in ast.walk(st):
                if isinstance(n, ast.Name) and n.id == name \
                        and isinstance(n.ctx, ast.Load):
                    return True
        return False

    root = Path(__file__).resolve().parent.parent
    offenders = []
    for f in sorted(list((root / "scripts").glob("*.py"))
                    + list((root / "src").rglob("*.py"))):
        try:
            tree = ast.parse(f.read_text())
        except SyntaxError:
            continue
        for outer in ast.walk(tree):
            if not isinstance(outer, ast.For):
                continue
            ot = targets(outer)
            for inner in ast.walk(outer):
                if inner is outer or not isinstance(inner, ast.For):
                    continue
                for name in ot & targets(inner):
                    if read_after(outer, name, inner):
                        offenders.append(
                            f"{f.relative_to(root)}:{inner.lineno} rebinds "
                            f"'{name}' from the loop at line {outer.lineno}")
    assert not offenders, "loop-variable shadowing:\n  " + "\n  ".join(offenders)
