import pandas as pd
import pytest
from enzsolv.composition_ablation import assign_saved_folds,composition_model_spec


def test_composition_ablation_omits_family_but_retains_conditions_metadata_and_composition():
    numeric,categorical=composition_model_spec(pd.DataFrame(columns=['aa_A','aa_C','identity30','target_log1p']))
    assert 'aa_A' in numeric and 'aa_C' in numeric and 'solvent_volume' in numeric
    assert 'meta_control_water' in numeric
    assert categorical==['solvent_name']
    assert 'identity30' not in numeric+categorical and 'target_log1p' not in numeric+categorical


def test_saved_folds_are_matched_by_id_and_require_full_unique_coverage():
    frame=pd.DataFrame({'measurement_id':[2,1,3],'sequence_id':['b','a','c'],'doi':['q','p','r']})
    assignments=pd.DataFrame({'measurement_id':[1,2,3],'fold':[0,1,2]})
    result=assign_saved_folds(frame,assignments)
    assert result.fold.tolist()==[1,0,2]
    assert 'fold' not in frame
    with pytest.raises(ValueError,match='coverage'):
        assign_saved_folds(frame,assignments.iloc[:2])
    with pytest.raises(ValueError,match='Duplicate'):
        assign_saved_folds(frame,pd.concat([assignments,assignments.iloc[:1]]))


def test_saved_folds_reject_publication_or_sequence_crossing():
    frame=pd.DataFrame({'measurement_id':[1,2],'sequence_id':['a','b'],'doi':['p','p']})
    assignments=pd.DataFrame({'measurement_id':[1,2],'fold':[0,1]})
    with pytest.raises(ValueError,match='publication'):
        assign_saved_folds(frame,assignments)
