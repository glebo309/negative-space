"""Exact comparator-matched ablation requested by final independent review."""
from .run_benchmark import CONDITIONS
from .metadata_diagnostics import META_FEATURES


def composition_model_spec(frame):
    return CONDITIONS+META_FEATURES+sorted(c for c in frame if c.startswith('aa_')),['solvent_name']


def assign_saved_folds(frame,assignments):
    if frame.measurement_id.duplicated().any() or assignments.measurement_id.duplicated().any():
        raise ValueError('Duplicate measurement identifier')
    if set(frame.measurement_id)!=set(assignments.measurement_id):
        raise ValueError('Saved fold coverage differs from training cohort')
    result=frame.copy()
    result['fold']=result.measurement_id.map(assignments.set_index('measurement_id').fold)
    if result.fold.isna().any():
        raise ValueError('Missing fold')
    if result.groupby('sequence_id').fold.nunique().max()!=1:
        raise ValueError('A sequence crosses folds')
    if result.groupby('doi').fold.nunique().max()!=1:
        raise ValueError('A publication crosses folds')
    return result
