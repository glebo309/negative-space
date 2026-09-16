import numpy as np
import pandas as pd

from enzsolv.features import AMINO_ACIDS, add_amino_acid_composition, grouped_folds


def test_add_amino_acid_composition_returns_normalized_frequencies():
    frame = pd.DataFrame({"sequence": ["AACC", "GGGG"]})

    featured = add_amino_acid_composition(frame)

    assert featured.loc[0, "aa_A"] == 0.5
    assert featured.loc[0, "aa_C"] == 0.5
    assert featured.loc[1, "aa_G"] == 1.0
    assert np.allclose(featured[[f"aa_{aa}" for aa in AMINO_ACIDS]].sum(axis=1), 1)


def test_grouped_folds_keep_sequences_out_of_training():
    frame = pd.DataFrame(
        {
            "sequence": ["AAAA", "AAAA", "BBBB", "BBBB", "CCCC", "CCCC"],
            "value": range(6),
        }
    )

    folds = list(grouped_folds(frame, n_splits=3))

    assert len(folds) == 3
    for train_index, test_index in folds:
        train_sequences = set(frame.iloc[train_index]["sequence"])
        test_sequences = set(frame.iloc[test_index]["sequence"])
        assert train_sequences.isdisjoint(test_sequences)
