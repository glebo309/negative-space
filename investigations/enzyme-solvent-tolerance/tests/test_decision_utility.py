import pandas as pd
import pytest
from enzsolv.decision_utility import utility_scores


def test_utility_averages_series_then_sequences_and_handles_prediction_ties():
    frame=pd.DataFrame({'sequence_id':['a']*6,'doi':['p']*3+['q']*3,
        'cluster':['g']*6,'solvent_name':['X','Y','Z']*2,
        'measured_value':[0.,50.,100.]*2,
        'candidate':[0.,1.,2.,0.,0.,0.],'reference':[2.,1.,0.,0.,1.,2.]})
    scores=utility_scores(frame,'candidate','reference')
    assert len(scores)==1
    s=scores.iloc[0]
    assert s.n_series==2
    assert s.candidate_pairwise==.75  # perfect then all ties
    assert s.reference_pairwise==.5  # reversed then perfect
    assert s.candidate_top_rank==.75  # best then random tie across all three
    assert s.reference_top_rank==.5
    assert s.candidate_best_hit==pytest.approx(2/3)
    assert s.reference_best_hit==.5


def test_utility_ignores_observed_tie_pairs_and_averages_duplicate_solvents():
    frame=pd.DataFrame({'sequence_id':['a']*4,'cluster':['g']*4,
        'solvent_name':['X','X','Y','Z'],'measured_value':[0.,20.,10.,20.],
        'candidate':[0.,0.,0.,1.],'reference':[0.,0.,1.,0.]})
    s=utility_scores(frame,'candidate','reference').iloc[0]
    assert s.candidate_pairwise==1.
    assert s.reference_pairwise==.25  # X-Z tied predictions, Y-Z reversed
    assert s.candidate_top_rank==1.
    assert s.reference_top_rank==.25


def test_utility_drops_panels_with_too_few_solvents_or_constant_response():
    frame=pd.DataFrame({'sequence_id':['a']*3+['b']*2,'cluster':['g']*5,
        'solvent_name':['X','Y','Z','X','Y'],'measured_value':[1.,1.,1.,0.,1.],
        'candidate':[1.]*5,'reference':[1.]*5})
    assert utility_scores(frame,'candidate','reference').empty
