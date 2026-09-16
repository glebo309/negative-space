"""Decision-facing ranking diagnostics from saved out-of-fold predictions."""
import numpy as np
import pandas as pd


def utility_scores(frame,candidate,reference):
    keys=['sequence_id']+[c for c in ['doi','solvent_volume','procedure','temperature','ph',
        'assay_solution','substrates','substrate_concentrations','cofactor','shaking','comments'] if c in frame]
    if frame.groupby('sequence_id').cluster.nunique().max()>1:
        raise ValueError('A sequence has multiple uncertainty groups')
    records=[]
    for _,group in frame.groupby(keys,dropna=False):
        panel=group.groupby('solvent_name')[['measured_value',candidate,reference]].mean()
        if len(panel)<3 or panel.measured_value.nunique()<2:
            continue
        observed=panel.measured_value.to_numpy()
        percentile=(panel.measured_value.rank(method='average').to_numpy()-1)/(len(panel)-1)
        i,j=np.triu_indices(len(panel),1)
        direction=np.sign(observed[i]-observed[j])
        valid=direction!=0
        record=dict(sequence_id=group.sequence_id.iloc[0],cluster=group.cluster.iloc[0])
        for label,column in [('candidate',candidate),('reference',reference)]:
            predicted=panel[column].to_numpy()
            if not np.isfinite(predicted).all():
                raise ValueError('Non-finite predictions')
            comparison=np.sign(predicted[i]-predicted[j])
            credit=np.where(comparison[valid]==0,.5,(comparison[valid]==direction[valid]).astype(float))
            selected=predicted==predicted.max()
            record[label+'_pairwise']=float(credit.mean())
            record[label+'_top_rank']=float(percentile[selected].mean())
            record[label+'_best_hit']=float((observed[selected]==observed.max()).mean())
        records.append(record)
    metrics=[label+'_'+metric for label in ['candidate','reference'] for metric in ['pairwise','top_rank','best_hit']]
    if not records:
        return pd.DataFrame(columns=['sequence_id','cluster','n_series']+metrics)
    scored=pd.DataFrame(records)
    result=scored.groupby(['sequence_id','cluster'],as_index=False)[metrics].mean()
    counts=scored.groupby(['sequence_id','cluster']).size().rename('n_series').reset_index()
    return result.merge(counts,on=['sequence_id','cluster'],validate='one_to_one')
