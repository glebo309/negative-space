"""Auditable literal metadata flags, not a claim of source-level curation."""

META_FEATURES=['meta_control_nonincubated','meta_control_water','meta_control_unknown',
               'meta_assay_aqueous','meta_assay_organic','meta_mode_uncertain',
               'meta_organic_reference']


def annotate_metadata(frame):
    result=frame.copy()
    comments=result.comments.fillna('').str.lower()
    final=comments.str.split('|').str[-1]
    procedure=result.procedure.fillna('').str.lower()
    result['meta_control_nonincubated']=comments.str.contains('non-incubated control',regex=False).astype(int)
    result['meta_control_water']=comments.str.contains('water-incubated control',regex=False).astype(int)
    result['meta_control_unknown']=(result.meta_control_nonincubated.eq(0)&result.meta_control_water.eq(0)).astype(int)
    result['meta_assay_aqueous']=procedure.str.contains('in aqueous phase',regex=False).astype(int)
    result['meta_assay_organic']=procedure.str.contains('in the presence of organic solvent',regex=False).astype(int)
    result['meta_mode_uncertain']=(final.str.contains('uncertain',regex=False)&
                                   final.str.contains('assay|activity|control',regex=True)).astype(int)
    result['meta_organic_reference']=comments.str.contains(r'control[^|]*presence of organic solvent',regex=True).astype(int)
    return result


def cohort_masks(frame):
    qualified=frame.meta_mode_uncertain.eq(0)&frame.meta_organic_reference.eq(0)
    return {'all':frame.index.to_series().notna(), 'mode_qualified':qualified,
            'aqueous_explicit':qualified&frame.meta_assay_aqueous.eq(1)}
