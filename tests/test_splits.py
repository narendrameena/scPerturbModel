import numpy as np
import pandas as pd

from perturbmodel.evaluation import held_out_condition_triples, split_column


def toy_df():
    rng = np.random.default_rng(42)
    lines = [f"L{i}" for i in range(6)]
    drugs = [f"D{i}" for i in range(20)] + ["DMSO_TF"]
    rows = [(l, d, c) for l in lines for d in drugs
            for c in ([0.05, 0.5, 5.0] if d != "DMSO_TF" else [0.0])
            if rng.random() > 0.1]
    return pd.DataFrame(rows, columns=["cell_line_id", "drug", "conc"])


def test_split_is_deterministic_and_order_invariant():
    df = toy_df()
    shuffled = df.sample(frac=1.0, random_state=7)
    assert held_out_condition_triples(df, seed=0) == \
           held_out_condition_triples(shuffled, seed=0)
    assert held_out_condition_triples(df, seed=0) != \
           held_out_condition_triples(df, seed=1)


def test_every_test_perturbation_is_seen_in_training():
    df = toy_df()
    test = held_out_condition_triples(df, seed=0, test_frac=0.3)
    assert test, "split produced no test conditions"
    train = {tuple(r) for r in df.itertuples(index=False)} - test
    train_pairs = {(d, c) for _, d, c in train}
    assert all((d, c) in train_pairs for _, d, c in test)
    assert all(d != "DMSO_TF" for _, d, c in test)


def test_split_column_labels():
    df = toy_df()
    test = held_out_condition_triples(df, seed=0)
    col = split_column(df, test, seed=0, val_frac=0.2)
    assert set(col.unique()) <= {"train", "val", "test"}
    is_test = np.array([tuple(r) in test for r in df.itertuples(index=False)])
    assert (col[is_test] == "test").all()
    assert (col[~is_test] != "test").all()
