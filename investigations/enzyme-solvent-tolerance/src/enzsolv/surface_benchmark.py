"""Matched-cohort comparisons for protein surface descriptors."""
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from .benchmark import paired_cluster_ci, rank_scores


def align_cohort(rows, surfaces):
    if surfaces.sequence_id.duplicated().any():
        raise ValueError('Duplicate surface sequence identifiers')
    return rows.merge(surfaces,on='sequence_id',how='inner',sort=False,validate='many_to_one')


def feature_sets(frame):
    surface = sorted(c for c in frame if c.startswith('surface_'))
    composition = sorted(c for c in frame if c.startswith('aa_'))
    esm = sorted(c for c in frame if c.startswith('esm_'))
    return dict(conditions=[],composition=composition,surface=surface,
                esm650=esm,esm650_surface=esm+surface)


def model_pipeline(numeric, seed=11, trees=250):
    return Pipeline([
        ('prep',ColumnTransformer([
            ('num',SimpleImputer(strategy='median',add_indicator=True,
                                 keep_empty_features=True),numeric),
            ('cat',OneHotEncoder(handle_unknown='ignore',sparse_output=False),['solvent_name'])])),
        ('model',ExtraTreesRegressor(n_estimators=trees,min_samples_leaf=5,
                  max_features=1.,n_jobs=4,random_state=seed))])


def evaluate(rows, candidate, reference, matched=True):
    compared = rows.copy()
    # Assign simultaneously: one of the requested columns may itself be composition.
    compared[['baseline','composition']] = rows[[reference,candidate]].to_numpy()
    scores = rank_scores(compared,matched=matched)
    statistics = paired_cluster_ci(scores)
    statistics.update(candidate=candidate,reference=reference,
        candidate_mean=float(scores.composition.mean()),reference_mean=float(scores.baseline.mean()),
        n_series=int(scores.n_series.sum()))
    return statistics,scores.rename(columns={'baseline':'reference_score','composition':'candidate_score'})
