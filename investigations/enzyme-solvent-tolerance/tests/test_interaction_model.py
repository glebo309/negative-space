import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge

from enzsolv.interaction_model import InteractionFeatures


def test_interactions_can_reverse_solvent_preference_between_proteins():
    frame=pd.DataFrame({'sequence_id':['p','p','q','q'], 'solvent_name':['A','B','A','B'],
                        'condition':[1.,1.,1.,1.], 'protein_feature':[-1.,-1.,1.,1.]})
    target=np.array([0.,1.,1.,0.])
    encoder=InteractionFeatures(['condition'],['protein_feature'])
    x=encoder.fit_transform(frame)
    prediction=Ridge(alpha=.001).fit(x,target).predict(x)
    assert prediction[0]<prediction[1] and prediction[2]>prediction[3]


def test_protein_scaling_uses_unique_training_sequences_and_unknown_solvents_are_safe():
    frame=pd.DataFrame({'sequence_id':['p','p','p','q'], 'solvent_name':['A']*4,
                        'condition':[1.,2.,3.,4.], 'protein_feature':[0.,0.,0.,10.]})
    encoder=InteractionFeatures(['condition'],['protein_feature'])
    fitted=encoder.fit_transform(frame)
    assert encoder.protein_scaler_.mean_[0]==5.
    assert encoder.context_.named_steps['scale'].mean_[0]==2.5
    unseen=pd.DataFrame({'sequence_id':['x'], 'solvent_name':['unseen'],
                        'condition':[np.nan], 'protein_feature':[1000.]})
    transformed=encoder.transform(unseen)
    assert transformed.shape[1]==fitted.shape[1]
    assert np.isfinite(transformed).all()
    assert encoder.protein_scaler_.mean_[0]==5.


def test_conditions_only_has_solvent_condition_interactions_and_no_protein_terms():
    frame=pd.DataFrame({'sequence_id':['p','q'], 'solvent_name':['A','B'], 'condition':[1.,2.]})
    encoder=InteractionFeatures(['condition'],[])
    transformed=encoder.fit_transform(frame)
    assert transformed.shape==(2,5)  # 1 condition + 2 solvent indicators + 2 products
