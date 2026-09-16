"""Study-independent easier tasks with explicit support and uncertainty units."""
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from .benchmark import connected_groups
from .run_benchmark import CONDITIONS


def uncertainty_groups(rows, unit):
    nodes=sorted(rows[unit].unique())
    edges=[]
    for _, group in rows.groupby('doi'):
        units=group[unit].unique()
        edges.extend((units[0],value) for value in units[1:])
    return connected_groups(nodes,edges)


def sequence_publication_folds(rows, n_splits=5, random_state=None):
    assigned=rows.copy()
    mapping=uncertainty_groups(rows,'sequence_id')
    assigned['split_group']=assigned.sequence_id.map(mapping)
    assigned['fold']=-1
    splitter=GroupKFold(n_splits=n_splits,shuffle=random_state is not None,random_state=random_state)
    for fold,(_,test) in enumerate(splitter.split(assigned,groups=assigned.split_group)):
        assigned.iloc[test,assigned.columns.get_loc('fold')]=fold
    return assigned


def familiar_support(rows, original, minimum_overall=5, minimum_train=3):
    overall=original.groupby('identity30').agg(sequences=('sequence_id','nunique'),papers=('doi','nunique'))
    qualified=set(overall.index[(overall.sequences>=minimum_overall)&(overall.papers>=minimum_overall)])
    records=[]
    for fold in sorted(rows.fold.unique()):
        train=rows.loc[rows.fold.ne(fold)]
        counts=train.groupby('identity30').agg(training_sequences=('sequence_id','nunique'),training_papers=('doi','nunique'))
        for family in rows.loc[rows.fold.eq(fold),'identity30'].unique():
            nseq,npub=(counts.loc[family].tolist() if family in counts.index else [0,0])
            records.append(dict(fold=fold,identity30=family,training_sequences=int(nseq),training_papers=int(npub),
                eligible=family in qualified and nseq>=minimum_train and npub>=minimum_train))
    return pd.DataFrame(records)


def known_publication_splits(rows):
    output=[]
    repeated=set(rows.groupby('sequence_id').doi.nunique().loc[lambda x:x>=2].index)
    for publication in sorted(rows.loc[rows.sequence_id.isin(repeated),'doi'].unique()):
        train=rows.doi.ne(publication)
        known=set(rows.loc[train,'sequence_id'])
        test=rows.doi.eq(publication)&rows.sequence_id.isin(known)
        if test.any():
            output.append(dict(publication=publication,train=train,test=test))
    return output


def assert_split(train, test, known=False):
    if set(train.doi)&set(test.doi):
        raise ValueError('A publication crosses train/test')
    trainseq,testseq=set(train.sequence_id),set(test.sequence_id)
    if known and not testseq<=trainseq:
        raise ValueError('A test sequence is absent from training')
    if not known and trainseq&testseq:
        raise ValueError('A sequence crosses train/test')


def model_pipeline(numeric, categorical, seed=11, trees=250):
    return Pipeline([
        ('prep',ColumnTransformer([
            ('num',SimpleImputer(strategy='median',add_indicator=True,keep_empty_features=True),numeric),
            ('cat',OneHotEncoder(handle_unknown='ignore',sparse_output=False),categorical)])),
        ('model',ExtraTreesRegressor(n_estimators=trees,min_samples_leaf=5,max_features=1.,
                                    n_jobs=4,random_state=seed))])


def familiar_model_specs(sets,metadata=None,names=None):
    base=CONDITIONS+([] if metadata is None else metadata)
    models=dict(conditions=(base,['solvent_name']),family_conditions=(base,['solvent_name','identity30']))
    for name in ['composition','surface','esm650']:
        models[name]=(base+sets[name],['solvent_name','identity30'])
    return models if names is None else {name:models[name] for name in names}


def familiar_pairs(models):
    names=set(models)
    pairs=[]
    if {'family_conditions','conditions'}<=names:
        pairs.append(('family_conditions','conditions'))
    for reference in ['family_conditions','conditions']:
        if reference in names:
            pairs += [(name,reference) for name in ['composition','surface','esm650'] if name in names]
    if 'composition' in names:
        pairs += [(name,'composition') for name in ['surface','esm650'] if name in names]
    return pairs


def aggregate_repeated_scores(tables):
    rows=pd.concat(tables,ignore_index=True)
    return rows.groupby('sequence_id',as_index=False).agg(
        candidate_score=('candidate_score','mean'),reference_score=('reference_score','mean'),
        available_runs=('candidate_score','size'),mean_series_per_run=('n_series','mean'))
