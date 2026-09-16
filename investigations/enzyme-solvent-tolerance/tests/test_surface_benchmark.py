import numpy as np
import pandas as pd
import pytest

from enzsolv.surface_benchmark import align_cohort, feature_sets, model_pipeline, evaluate


def test_align_cohort_preserves_saved_rows_and_folds_and_excludes_missing_structures():
    rows = pd.DataFrame({'measurement_id':[10,11,12], 'sequence_id':['b','a','c'],
        'sequence':['AD','AC','DE'], 'fold':[2,0,1], 'cluster':['b','a','c'],
        'measured_value':[10.,20.,30.]})
    surfaces = pd.DataFrame({'sequence_id':['a','b'], 'surface_aa_A':[.2,.7]})
    aligned = align_cohort(rows,surfaces)
    assert aligned.measurement_id.tolist()==[10,11]
    assert aligned.fold.tolist()==[2,0]
    assert aligned.surface_aa_A.tolist()==[.7,.2]
    with pytest.raises(ValueError,match='Duplicate'):
        align_cohort(rows,pd.concat([surfaces,surfaces.iloc[:1]]))


def test_feature_sets_are_explicit_and_do_not_include_identifiers_or_outcomes():
    frame = pd.DataFrame(columns=['sequence_id','measured_value','fold','cluster',
        'target_log1p','mean_plddt','surface_aa_A','aa_A','esm_000','solvent_volume'])
    sets = feature_sets(frame)
    assert sets['surface']==['surface_aa_A']
    assert sets['composition']==['aa_A']
    assert sets['esm650_surface']==['esm_000','surface_aa_A']
    assert sets['conditions']==[]


def test_model_pipeline_fits_imputation_only_on_training_rows():
    train = pd.DataFrame({'condition':[1.,3.,np.nan], 'surface_aa_A':[.1,.3,.2],
                          'solvent_name':['A','A','B']})
    model = model_pipeline(['condition','surface_aa_A'],seed=11,trees=5)
    model.fit(train,np.array([1.,2.,3.]))
    assert np.allclose(model.named_steps['prep'].named_transformers_['num'].statistics_,[2.,.2])
    test = pd.DataFrame({'condition':[1000.], 'surface_aa_A':[.9], 'solvent_name':['new']})
    assert np.isfinite(model.predict(test)).all()
    assert model.named_steps['prep'].named_transformers_['num'].statistics_[0]==2.


def test_evaluate_returns_paired_comparisons_on_identical_matched_series():
    rows = pd.DataFrame({'sequence_id':['s']*3, 'cluster':['s']*3,
        'solvent_name':['A','B','C'], 'measured_value':[1.,2.,3.],
        'conditions':[3.,2.,1.], 'surface':[1.,2.,3.]})
    statistics,scores = evaluate(rows,'surface','conditions',matched=True)
    assert statistics['candidate_mean']==1.
    assert statistics['reference_mean']==-1.
    assert statistics['delta']==2.
    assert len(scores)==1
