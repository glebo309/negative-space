import pandas as pd
import pytest

from enzsolv.decision_diagnostics import repeat_groups, support_tables, planted_target, matched_rank_targets, optimistic_row_folds


def test_planted_signal_changes_solvent_preference_without_mutating_inputs():
    frame=pd.DataFrame({'protein_signal':[-1.,-1.,1.,1.], 'solvent_logp':[-1.,1.,-1.,1.]})
    original=frame.copy(deep=True)
    strong=planted_target(frame,strength=1.5,noise=0)
    null=planted_target(frame,strength=0,noise=0)
    assert strong[0]>strong[1] and strong[2]<strong[3]
    assert null[0]<null[1] and null[2]<null[3]
    pd.testing.assert_frame_equal(frame,original)


def test_repeats_require_distinct_papers_and_complete_conditions():
    frame=pd.DataFrame({'sequence_id':['a']*5,'doi':['p','p','q','r','s'],
        'solvent_name':['A']*5,'temperature':[20.,20.,20.,None,None],
        'measured_value':[10.,12.,30.,100.,1.]})
    repeats=repeat_groups(frame,['sequence_id','solvent_name','temperature'])
    assert len(repeats)==1
    assert repeats.iloc[0].papers==2
    assert repeats.iloc[0].range_percentage_points==19.  # paper p mean 11 versus q 30
    assert repeats.iloc[0].rows==3


def test_support_counts_proteins_and_publications_not_measurement_rows():
    frame=pd.DataFrame({'measurement_id':[1,2,3,4], 'sequence_id':['a','a','b','c'],
        'doi':['p','q','q','r'], 'enzyme_name':['A','A','B','C'],
        'solvent_name':['X','Y','X','X'], 'family':['f','f','f','g']})
    sequences,families=support_tables(frame,'family')
    assert sequences.set_index('sequence_id').loc['a','papers']==2
    assert families.set_index('family').loc['f','sequences']==2
    assert families.set_index('family').loc['f','papers']==2
    assert families.set_index('family').loc['f','rows']==3
    with pytest.raises(ValueError,match='measurement'):
        support_tables(pd.concat([frame,frame.iloc[:1]]),'family')


def test_rank_targets_average_duplicate_solvents_and_exclude_unscoreable_series():
    frame=pd.DataFrame({'sequence_id':['a']*6+['b']*3+['c']*2,
        'doi':['p']*4+['q']*2+['p']*5,
        'solvent_name':['X','X','Y','Z','X','Y','X','Y','Z','X','Y'],
        'measured_value':[0.,20.,20.,30.,0.,100.,1.,1.,1.,0.,1.]},
        index=range(10,21))
    original=frame.copy(deep=True)
    target=matched_rank_targets(frame)
    assert target.index.equals(frame.index)
    assert target.iloc[:4].tolist()==[0.,0.,.5,1.]
    assert target.iloc[4:].isna().all()
    pd.testing.assert_frame_equal(frame,original)


def test_rank_targets_use_average_ranks_for_ties():
    frame=pd.DataFrame({'sequence_id':['a']*3,'solvent_name':['X','Y','Z'],
                        'measured_value':[10.,10.,20.]})
    assert matched_rank_targets(frame).tolist()==[.25,.25,1.]


def test_optimistic_row_folds_cover_each_row_once_but_deliberately_allow_sequence_overlap():
    import numpy as np
    frame=pd.DataFrame({'sequence_id':['a']*10+['b']*10})
    fold=optimistic_row_folds(frame,n_splits=5,seed=11)
    assert len(fold)==20 and sorted(set(fold))==[0,1,2,3,4]
    np.testing.assert_array_equal(fold,optimistic_row_folds(frame,n_splits=5,seed=11))
    for k in set(fold):
        assert set(frame.loc[fold==k,'sequence_id'])&set(frame.loc[fold!=k,'sequence_id'])
