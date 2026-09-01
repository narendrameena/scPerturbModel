import numpy as np
import pandas as pd
import pytest

from perturbmodel.decompose import decompose


def synth(n_ctx=6, n_pert=20, n_rep=3, n_genes=200, interaction=1.0, seed=0):
    """Atlas with a known additive effect plus a tunable context x perturbation
    interaction, so the estimator can be checked against ground truth."""
    rng = np.random.default_rng(seed)
    shared = {p: rng.normal(0, 1, n_genes) for p in range(n_pert)}
    inter = {(c, p): rng.normal(0, interaction, n_genes)
             for c in range(n_ctx) for p in range(n_pert)}
    rows, meta = [], []
    for c in range(n_ctx):
        base = rng.normal(5, 1, n_genes)
        for r in range(n_rep):
            rows.append(base + rng.normal(0, 0.3, n_genes))
            meta.append({"ctx": f"c{c}", "pert": "CTRL", "rep": f"r{r}"})
            for p in range(n_pert):
                rows.append(base + shared[p] + inter[(c, p)]
                            + rng.normal(0, 0.3, n_genes))
                meta.append({"ctx": f"c{c}", "pert": f"p{p}", "rep": f"r{r}"})
    import anndata as ad
    obs = pd.DataFrame(meta)
    obs.index = obs.index.astype(str)
    return ad.AnnData(X=np.stack(rows).astype(np.float32), obs=obs)


def test_recovers_a_known_interaction():
    a = synth(interaction=1.0)
    r = decompose(a, context="ctx", perturbation="pert", control="CTRL",
                  replicate="rep", n_genes=200)
    assert r.n_contexts == 6 and r.n_perturbations == 20
    assert r.interaction > 0
    assert 0.2 < r.interaction_share < 0.8


def test_interaction_shrinks_when_absent():
    strong = decompose(synth(interaction=1.5, seed=1), context="ctx",
                       perturbation="pert", control="CTRL", replicate="rep",
                       n_genes=200)
    none = decompose(synth(interaction=0.0, seed=1), context="ctx",
                     perturbation="pert", control="CTRL", replicate="rep",
                     n_genes=200)
    assert none.interaction_share < strong.interaction_share
    assert none.interaction_share < 0.15


def test_residual_reproduces_across_replicates_not_perturbations():
    r = decompose(synth(interaction=1.0), context="ctx", perturbation="pert",
                  control="CTRL", replicate="rep", n_genes=200)
    rc = r.residual_correlations.set_index("comparison").median_r
    across_rep = rc.get("across replicates (same ctx+pert)")
    across_pert = rc.get("across perturbations (same ctx)")
    assert across_rep is not None and across_pert is not None
    assert across_rep > across_pert          # the defining signature


def test_warns_without_replicates():
    a = synth(n_rep=1)
    r = decompose(a, context="ctx", perturbation="pert", control="CTRL",
                  n_genes=200)
    assert any("replicate" in w for w in r.warnings)


def test_rejects_missing_control():
    with pytest.raises(ValueError):
        decompose(synth(), context="ctx", perturbation="pert",
                  control="NOT_A_CONTROL", replicate="rep", n_genes=200)


def test_detects_a_designed_replicate_batch():
    """A batch that repeats another's conditions should be flagged.

    This is the failure that cost this project its headline number: Tahoe's
    plate 14 duplicates plate 6, is withheld from training by convention, and
    was therefore dropped by our own default -- which removed every same-dose
    replicate and doubled the apparent interaction share.
    """
    from perturbmodel.decompose import find_replicate_batches
    a = synth(n_ctx=4, n_pert=6, n_rep=2)
    obs = a.obs.copy()
    # rep r0/r1 measure the same ctx x pert conditions, so they are replicates
    found = find_replicate_batches(obs, "ctx", "pert", "rep")
    assert len(found) >= 1
    assert found.jaccard.iloc[0] > 0.9


def test_no_false_replicate_when_batches_are_disjoint():
    import pandas as pd
    from perturbmodel.decompose import find_replicate_batches
    a = synth(n_ctx=4, n_pert=6, n_rep=2)
    obs = a.obs.copy()
    # give each replicate its own disjoint perturbation set
    obs.loc[obs.rep == "r0", "pert"] = obs.loc[obs.rep == "r0", "pert"] + "_A"
    obs.loc[obs.rep == "r1", "pert"] = obs.loc[obs.rep == "r1", "pert"] + "_B"
    assert len(find_replicate_batches(obs, "ctx", "pert", "rep")) == 0
