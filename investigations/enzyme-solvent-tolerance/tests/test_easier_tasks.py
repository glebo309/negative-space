import numpy as np
import pandas as pd
import pytest

from enzsolv.easier_tasks import (sequence_publication_folds, familiar_support,
    known_publication_splits, uncertainty_groups, assert_split, model_pipeline)


def toy():
    return pd.DataFrame([dict(measurement_id=i, sequence_id=f's{i}', doi=f'p{i}',
        identity30='family', solvent_name='ethanol', x=float(i)) for i in range(6)])


def test_publication_and_sequence_components_never_cross_folds():
    rows=toy()
    rows=pd.concat([rows, rows.iloc[[0]].assign(measurement_id=7,doi='p1')],ignore_index=True)
    assigned=sequence_publication_folds(rows,n_splits=3)
    assert assigned.groupby('sequence_id').fold.nunique().max()==1
    assert assigned.groupby('doi').fold.nunique().max()==1
    assert assigned.loc[assigned.sequence_id.isin(['s0','s1']),'fold'].nunique()==1


def test_family_support_counts_distinct_training_units_and_prequalification():
    rows=toy().assign(fold=[0,0,1,1,2,2])
    support=familiar_support(rows,rows)
    assert support.eligible.all()
    assert support.training_sequences.eq(4).all()
    duplicated=pd.concat([rows,rows.iloc[[0]]]*3,ignore_index=True)
    assert familiar_support(duplicated,rows).training_sequences.eq(4).all()
    sparse=rows.iloc[:4]
    assert not familiar_support(sparse,sparse).eligible.any()


def test_known_publication_splits_hold_entire_paper_and_require_seen_sequence():
    rows=toy()
    rows=pd.concat([rows,rows.iloc[[0]].assign(measurement_id=7,doi='p1')],ignore_index=True)
    splits=known_publication_splits(rows)
    assert {s['publication'] for s in splits}=={'p0','p1'}
    for s in splits:
        train=rows.loc[s['train']]
        test=rows.loc[s['test']]
        assert set(train.doi).isdisjoint(test.doi)
        assert set(test.sequence_id)<=set(train.sequence_id)


def test_uncertainty_connects_families_sharing_papers_without_changing_folds():
    rows=toy().assign(identity30=['a','b','c','c','d','d'],doi=['p','p','q','r','s','t'])
    groups=uncertainty_groups(rows,'identity30')
    assert groups['a']==groups['b']
    assert groups['a']!=groups['c']


def test_split_assertions_distinguish_known_and_unseen_sequences():
    train=toy().iloc[:3]
    test=train.iloc[[0]].assign(doi='new')
    assert_split(train,test,known=True)
    with pytest.raises(ValueError,match='sequence'):
        assert_split(train,test,known=False)
    with pytest.raises(ValueError,match='publication'):
        assert_split(train,train.iloc[[0]],known=True)


def test_model_categories_fit_only_training_and_predictions_finite():
    rows=toy()
    model=model_pipeline(['x'],['solvent_name','identity30'],trees=3)
    model.fit(rows,np.arange(len(rows)))
    pred=model.predict(rows.iloc[[0]].assign(identity30='unseen',solvent_name='new'))
    assert np.isfinite(pred).all()
    assert 'unseen' not in model.named_steps['prep'].named_transformers_['cat'].categories_[1]


def test_fit_prediction_cache_binds_targets_and_keeps_test_target_unused(tmp_path):
    from enzsolv.run_easier_tasks import fit_prediction
    train=toy().assign(target_log1p=np.arange(6)/10)
    test=train.iloc[:2].copy().assign(doi='new',target_log1p=999)
    cache=tmp_path/'fit.npz'
    first=fit_prediction(train,test,['x'],['solvent_name'],cache,trees=3)
    again=fit_prediction(train,test.assign(target_log1p=-999),['x'],['solvent_name'],cache,trees=3)
    assert np.array_equal(first,again)
    assert np.isfinite(first).all()
    with pytest.raises(ValueError,match='signature'):
        fit_prediction(train.assign(target_log1p=100),test,['x'],['solvent_name'],cache,trees=3)


def test_second_assignment_reproducible_distinct_and_still_study_independent():
    rows=pd.concat([toy().assign(measurement_id=lambda x:x.measurement_id+10*j,
        sequence_id=lambda x:x.sequence_id+str(j),doi=lambda x:x.doi+str(j)) for j in range(3)],ignore_index=True)
    first=sequence_publication_folds(rows)
    second=sequence_publication_folds(rows,random_state=47)
    again=sequence_publication_folds(rows,random_state=47)
    assert second.fold.equals(again.fold)
    assert not first.fold.equals(second.fold)
    assert second.groupby('doi').fold.nunique().max()==1
    assert second.groupby('sequence_id').fold.nunique().max()==1


def test_robustness_model_features_and_baseline_comparisons():
    from enzsolv.easier_tasks import familiar_model_specs, familiar_pairs
    sets={'composition':['aa_A'],'surface':['surface_x'],'esm650':['esm_x']}
    models=familiar_model_specs(sets,metadata=['meta_x'],names=['family_conditions','composition','surface'])
    assert set(models)=={'family_conditions','composition','surface'}
    assert all('meta_x' in numeric for numeric,_ in models.values())
    assert all(categorical==['solvent_name','identity30'] for _,categorical in models.values())
    primary=familiar_model_specs(sets)
    pairs=familiar_pairs(primary)
    assert ('composition','conditions') in pairs
    assert ('surface','conditions') in pairs
    assert ('composition','family_conditions') in pairs
    assert ('esm650','composition') in pairs


def test_repeated_scores_average_within_sequence_and_count_each_once():
    from enzsolv.easier_tasks import aggregate_repeated_scores
    first=pd.DataFrame({'sequence_id':['a','b'],'candidate_score':[0.,0.],
        'reference_score':[0.,0.],'n_series':[2,3]})
    second=pd.DataFrame({'sequence_id':['a'],'candidate_score':[1.],
        'reference_score':[0.],'n_series':[2]})
    aggregate=aggregate_repeated_scores([first,second])
    assert len(aggregate)==2
    assert aggregate.candidate_score.mean()==.25
    assert aggregate.set_index('sequence_id').loc['a','available_runs']==2
    assert aggregate.set_index('sequence_id').loc['b','available_runs']==1
