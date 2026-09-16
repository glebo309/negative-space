import pandas as pd
from enzsolv.metadata_diagnostics import annotate_metadata, cohort_masks


def test_metadata_flags_are_literal_and_uncertainty_is_mode_specific():
    frame=pd.DataFrame({
        'comments':[
            'Non-incubated control | Uncertainty about substrate, here taken from reference',
            'Water-incubated control | Uncertainty about assay in aqueous phase',
            'Non-incubated control in the presence of organic solvent | Clear',
            None,
            'Uncertain sequence | Activity measured in aqueous phase'],
        'procedure':['Incubation in aqueous phase']*5})
    original=frame.copy(deep=True)
    result=annotate_metadata(frame)
    assert result.meta_mode_uncertain.tolist()==[0,1,0,0,0]
    assert result.meta_organic_reference.tolist()==[0,0,1,0,0]
    assert result.meta_control_water.tolist()==[0,1,0,0,0]
    assert result.meta_control_unknown.tolist()==[0,0,0,1,1]
    masks=cohort_masks(result)
    assert masks['mode_qualified'].tolist()==[True,False,False,True,True]
    assert masks['aqueous_explicit'].tolist()==[True,False,False,True,True]
    pd.testing.assert_frame_equal(frame,original)


def test_strict_filter_does_not_infer_aqueous_from_missing_or_organic_procedure():
    frame=pd.DataFrame({'comments':['Water-incubated control']*3,
                        'procedure':[None,'In the presence of organic solvent','In aqueous phase']})
    result=annotate_metadata(frame)
    assert cohort_masks(result)['aqueous_explicit'].tolist()==[False,False,True]
    assert result.meta_assay_organic.tolist()==[0,1,0]
