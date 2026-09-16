"""Support and controls for diagnosing the solvent-prediction benchmark."""
import numpy as np
import pandas as pd


def support_tables(frame,family):
    if frame.measurement_id.duplicated().any():
        raise ValueError('Duplicate measurement identifiers')
    sequences=frame.groupby('sequence_id',as_index=False).agg(
        rows=('measurement_id','size'),papers=('doi','nunique'),
        solvents=('solvent_name','nunique'),name=('enzyme_name','first'))
    families=frame.groupby(family,as_index=False).agg(
        rows=('measurement_id','size'),papers=('doi','nunique'),
        sequences=('sequence_id','nunique'),solvents=('solvent_name','nunique'),
        names=('enzyme_name',lambda x:'; '.join(x.unique()[:3])))
    return sequences,families


def repeat_groups(frame,keys):
    complete=frame.dropna(subset=keys+['doi','measured_value'])
    records=[]
    for values,group in complete.groupby(keys,dropna=False):
        paper_means=group.groupby('doi').measured_value.mean()
        if len(paper_means)<2:
            continue
        if not isinstance(values,tuple):
            values=(values,)
        records.append(dict(zip(keys,values),rows=len(group),papers=len(paper_means),
            range_percentage_points=float(paper_means.max()-paper_means.min()),
            minimum_paper_mean=float(paper_means.min()),maximum_paper_mean=float(paper_means.max()),
            dois='; '.join(str(x) for x in paper_means.index)))
    return pd.DataFrame(records,columns=keys+['rows','papers','range_percentage_points',
        'minimum_paper_mean','maximum_paper_mean','dois'])


def planted_target(frame,strength=1.5,noise=.1,seed=11):
    solvent=np.tanh(frame.solvent_logp.fillna(0).to_numpy()/2)
    protein=np.clip(frame.protein_signal.to_numpy(),-2,2)
    rng=np.random.default_rng(seed)
    return 5+.3*solvent+strength*protein*solvent+rng.normal(0,noise,len(frame))


def matched_rank_targets(frame):
    """Within-experiment solvent midranks; missing when the real metric is ineligible."""
    keys=['sequence_id']+[c for c in ['doi','solvent_volume','procedure','temperature','ph',
        'assay_solution','substrates','substrate_concentrations','cofactor','shaking','comments'] if c in frame]
    target=pd.Series(np.nan,index=frame.index,dtype=float)
    for _,group in frame.groupby(keys,dropna=False):
        observed=group.groupby('solvent_name').measured_value.mean()
        if len(observed)<3 or observed.nunique()<2:
            continue
        ranks=(observed.rank(method='average')-1)/(len(observed)-1)
        target.loc[group.index]=group.solvent_name.map(ranks)
    return target


def optimistic_row_folds(frame,n_splits=5,seed=11):
    """Deliberately leaky across proteins/studies: diagnostic only, never transfer CV."""
    from sklearn.model_selection import KFold
    folds=np.full(len(frame),-1,dtype=int)
    for fold,(_,test) in enumerate(KFold(n_splits=n_splits,shuffle=True,random_state=seed).split(frame)):
        folds[test]=fold
    return folds
