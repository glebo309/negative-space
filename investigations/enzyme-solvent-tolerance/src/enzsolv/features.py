import pandas as pd
from sklearn.model_selection import GroupKFold


AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")


def add_amino_acid_composition(frame: pd.DataFrame) -> pd.DataFrame:
    featured = frame.copy()
    sequences = featured["sequence"].str.upper()
    lengths = sequences.str.len()
    for amino_acid in AMINO_ACIDS:
        featured[f"aa_{amino_acid}"] = sequences.str.count(amino_acid) / lengths
    return featured


def grouped_folds(frame: pd.DataFrame, n_splits: int = 5):
    splitter = GroupKFold(n_splits=n_splits)
    return splitter.split(frame, groups=frame["sequence"])
