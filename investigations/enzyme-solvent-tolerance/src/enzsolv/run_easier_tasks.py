"""Run the pre-specified easier-task diagnostics without held-out tuning."""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import sklearn

from .benchmark import prepare_cohort, rank_scores, paired_cluster_ci
from .data import load_measurements
from .easier_tasks import (sequence_publication_folds, familiar_support,
    known_publication_splits, uncertainty_groups, assert_split, model_pipeline,
    familiar_model_specs, familiar_pairs, aggregate_repeated_scores)
from .embeddings import join_embeddings
from .features import add_amino_acid_composition
from .run_benchmark import CONDITIONS
from .surface_benchmark import align_cohort, feature_sets, evaluate

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/decision_diagnostics/easier_tasks'


def fit_prediction(train,test,numeric,categorical,cache,trees=250):
    features=numeric+categorical
    trainkeys=['measurement_id','sequence_id','doi','target_log1p']+features
    testkeys=['measurement_id','sequence_id','doi']+features
    trainkeys=list(dict.fromkeys(trainkeys))
    testkeys=list(dict.fromkeys(testkeys))
    signature=hashlib.sha256(
        pd.util.hash_pandas_object(train[trainkeys],index=False).values.tobytes()+
        pd.util.hash_pandas_object(test[testkeys],index=False).values.tobytes()+
        json.dumps(dict(numeric=numeric,categorical=categorical,trees=trees,seed=11,sklearn=sklearn.__version__)).encode()+
        inspect.getsource(model_pipeline).encode()+inspect.getsource(fit_prediction).encode()).hexdigest()
    if cache.exists():
        saved=np.load(cache,allow_pickle=False)
        if saved['signature'].item()!=signature:
            raise ValueError('Cache signature mismatch: '+str(cache))
        return saved['predictions']
    model=model_pipeline(numeric,categorical,trees=trees)
    model.fit(train,train.target_log1p)
    predictions=np.expm1(model.predict(test))
    assert np.isfinite(predictions).all()
    cache.parent.mkdir(parents=True,exist_ok=True)
    temporary=cache.with_suffix('.tmp.npz')
    np.savez_compressed(temporary,signature=np.array(signature),predictions=predictions)
    temporary.replace(cache)
    print('FIT',cache.relative_to(OUT) if OUT in cache.parents else cache.name,'complete',flush=True)
    return predictions


def slim(rows):
    return rows[[c for c in rows if not c.startswith(('esm_','aa_','surface_','confident_surface_')) and c!='sequence']].copy()


def support_summary(rows,unit):
    dummy=rows.copy().assign(baseline=0.,composition=0.,cluster=rows[unit])
    scores=rank_scores(dummy,matched=True)
    return dict(rows=len(rows),sequences=int(rows.sequence_id.nunique()),publications=int(rows.doi.nunique()),
        uncertainty_units=int(rows[unit].nunique()),scoreable_sequences=len(scores),series=int(scores.n_series.sum()))


def comparisons(predictions,pairs,folder,family=False):
    result={}
    seqfamily=predictions[['sequence_id','identity30']].drop_duplicates()
    for candidate,reference in pairs:
        label=candidate+'_vs_'+reference
        stats,scores=evaluate(predictions,candidate,reference,matched=True)
        scores=scores.merge(seqfamily,on='sequence_id',validate='one_to_one')
        scores['delta']=scores.candidate_score-scores.reference_score
        scores.to_csv(folder/f'{label}_scores.csv',index=False)
        if family:
            table=scores.groupby('identity30').agg(sequences=('sequence_id','size'),series=('n_series','sum'),
                candidate=('candidate_score','mean'),reference=('reference_score','mean'),delta=('delta','mean'))
            table.to_csv(folder/f'{label}_per_family.csv')
            influence=[]
            for omitted in sorted(scores.identity30.unique()):
                kept=scores.loc[scores.identity30.ne(omitted)]
                influence.append(dict(omitted_family=omitted,sequences=len(kept),
                    delta=float(kept.delta.mean()),candidate=float(kept.candidate_score.mean()),reference=float(kept.reference_score.mean())))
            pd.DataFrame(influence).to_csv(folder/f'{label}_leave_one_family_out_scores_only.csv',index=False)
        result[label]=stats
        print('RESULT',folder.name,label,json.dumps(stats),flush=True)
    return result


def familiar_task(rows,full,variant='familiar_family',metadata=None,names=None):
    folder=OUT/variant
    folder.mkdir(parents=True,exist_ok=True)
    support=familiar_support(rows,full)
    support.to_csv(folder/'train_test_support.csv',index=False)
    assigned=rows.merge(support,on=['fold','identity30'],validate='many_to_one')
    assigned[['measurement_id','sequence_id','doi','identity30','split_group','fold',
        'training_sequences','training_papers','eligible']].to_csv(folder/'fold_assignments.csv',index=False)
    eligible=assigned.loc[assigned.eligible].copy()
    groups=uncertainty_groups(eligible,'identity30')
    eligible['cluster']=eligible.identity30.map(groups)
    summary=support_summary(eligible,'cluster')
    summary['families']=int(eligible.identity30.nunique())
    (folder/'support_before_scores.json').write_text(json.dumps(summary,indent=2))
    print('SUPPORT familiar',json.dumps(summary),flush=True)
    sets=feature_sets(rows)
    models=familiar_model_specs(sets,metadata=metadata,names=names)
    predictions=slim(eligible)
    for name,(numeric,categorical) in models.items():
        predictions[name]=np.nan
        for fold in sorted(rows.fold.unique()):
            train=rows.loc[rows.fold.ne(fold)]
            test=eligible.loc[eligible.fold.eq(fold)]
            assert_split(train,test,known=False)
            assert set(test.identity30)<=set(train.identity30)
            if not len(test):
                continue
            values=fit_prediction(train,test,numeric,categorical,folder/f'{name}_fold{fold}.npz')
            predictions.loc[test.index,name]=values
        assert np.isfinite(predictions[name]).all()
        predictions.to_csv(folder/'predictions.csv',index=False)
    pairs=familiar_pairs(models)
    summary['comparisons']=comparisons(predictions,pairs,folder,family=True)
    return summary


def known_task(rows):
    folder=OUT/'known_enzyme'
    folder.mkdir(parents=True,exist_ok=True)
    splits=known_publication_splits(rows)
    selected=rows.loc[np.logical_or.reduce([s['test'] for s in splits])].copy()
    mapping=uncertainty_groups(selected,'sequence_id')
    selected['cluster']=selected.sequence_id.map(mapping)
    summary=support_summary(selected,'cluster')
    (folder/'support_before_scores.json').write_text(json.dumps(summary,indent=2))
    print('SUPPORT known',json.dumps(summary),flush=True)
    support=[]
    assignments=[]
    for fold,split in enumerate(splits):
        train=rows.loc[split['train']]
        test=rows.loc[split['test']]
        for sid,group in test.groupby('sequence_id'):
            training=train.loc[train.sequence_id.eq(sid)]
            support.append(dict(fold=fold,publication=split['publication'],sequence_id=sid,test_rows=len(group),
                training_sequence_rows=len(training),training_sequence_publications=int(training.doi.nunique()),
                training_sequence_solvents=int(training.solvent_name.nunique())))
        assignment=rows[['measurement_id','sequence_id','doi']].copy()
        assignment['fold']=fold
        assignment['role']=np.where(split['train'],'train',np.where(split['test'],'test_scored_candidate','test_unseen_excluded'))
        assignments.append(assignment)
    pd.DataFrame(support).to_csv(folder/'train_test_support.csv',index=False)
    pd.concat(assignments,ignore_index=True).to_csv(folder/'fold_assignments.csv',index=False)
    sets=feature_sets(rows)
    models=dict(conditions=(CONDITIONS,['solvent_name']),
        sequence_identity=(CONDITIONS,['solvent_name','sequence_id']),
        surface=(CONDITIONS+sets['surface'],['solvent_name']))
    predictions=slim(selected)
    for name,(numeric,categorical) in models.items():
        predictions[name]=np.nan
        for fold,split in enumerate(splits):
            train=rows.loc[split['train']]
            test=rows.loc[split['test']]
            assert_split(train,test,known=True)
            values=fit_prediction(train,test,numeric,categorical,folder/f'{name}_fold{fold}.npz')
            predictions.loc[test.index,name]=values
        assert np.isfinite(predictions[name]).all()
        predictions.to_csv(folder/'predictions.csv',index=False)
    summary['comparisons']=comparisons(predictions,[('sequence_identity','conditions'),
        ('surface','conditions'),('surface','sequence_identity')],folder)
    summary['leave_publication_fits_per_model']=len(splits)
    return summary


def run():
    started=time.monotonic()
    OUT.mkdir(parents=True,exist_ok=True)
    raw_path=ROOT/'data/raw/measurements.csv'
    expected='400823131c996a3899d95cdb0c82e37a6aa5833d32f0daa23e15a696045b3a06'
    assert hashlib.sha256(raw_path.read_bytes()).hexdigest()==expected
    full,_=prepare_cohort(load_measurements(raw_path))
    clusters=pd.read_csv(ROOT/'reports/family_benchmark/clusters.csv',usecols=['sequence_id','identity30'])
    full=full.merge(clusters,on='sequence_id',validate='many_to_one')
    full=sequence_publication_folds(full)
    full[['measurement_id','sequence_id','doi','identity30','split_group','fold']].to_csv(OUT/'full_cohort_folds.csv',index=False)
    rows=add_amino_acid_composition(full)
    saved=np.load(ROOT/'reports/esm650m/embeddings.npz',allow_pickle=False)
    seq=rows[['sequence_id','sequence']].drop_duplicates().sort_values('sequence_id')
    assert np.array_equal(saved['ids'],seq.sequence_id.to_numpy())
    assert np.array_equal(saved['hashes'],[hashlib.sha256(s.encode()).hexdigest() for s in seq.sequence])
    assert saved['revision'].item()=='08e4846e537177426273712802403f7ba8261b6c'
    rows=join_embeddings(rows,saved['ids'].tolist(),saved['vectors'])
    surfaces=pd.read_csv(ROOT/'reports/surface/features.csv')
    rows=align_cohort(rows,surfaces[['sequence_id']+[c for c in surfaces if c.startswith('surface_')]])
    assert rows.measurement_id.is_unique
    plan=dict(raw_sha256=expected,trees=250,min_samples_leaf=5,seed=11,target='log1p(percent)',
        family_qualification='at least 5 sequences and 5 publications in original cohort; 3 training sequences and 3 training papers per test family',
        family_folds='five-fold GroupKFold on full-cohort sequence/publication connected components before structure filter',
        known_folds='leave entire publication out, predict only exact sequences seen in another training publication',
        family_uncertainty='family/publication connected evaluation components, conditional paired bootstrap 5000 draws',
        known_uncertainty='sequence/publication connected evaluation components, conditional paired bootstrap 5000 draws',
        preprocess='training-only imputation and categorical encoding',
        inference='Exploratory, limited independent support; no held-out model tuning; identity tree permits interactions',
        score='Matched solvent Spearman, averaged by sequence',
        influence='leave-one-family-out score summaries only, no refitting')
    (OUT/'analysis_plan.json').write_text(json.dumps(plan,indent=2))
    result=dict(familiar_family=familiar_task(rows,full))
    (OUT/'results.json').write_text(json.dumps(result,indent=2))
    result['known_enzyme']=known_task(rows)
    result['checks']=dict(raw_sha256_after=hashlib.sha256(raw_path.read_bytes()).hexdigest(),
        unique_measurement_ids=bool(rows.measurement_id.is_unique),elapsed_seconds=time.monotonic()-started)
    assert result['checks']['raw_sha256_after']==expected
    (OUT/'results.json').write_text(json.dumps(result,indent=2))
    print('COMPLETE',json.dumps(result),flush=True)


def run_robustness():
    from .metadata_diagnostics import annotate_metadata, META_FEATURES
    expected=json.loads((OUT/'analysis_plan.json').read_text())['raw_sha256']
    raw_path=ROOT/'data/raw/measurements.csv'
    assert hashlib.sha256(raw_path.read_bytes()).hexdigest()==expected
    assert (OUT/'ROBUSTNESS_AMENDMENT.md').exists()
    primary_path=OUT/'familiar_family/predictions.csv'
    primary_hash=hashlib.sha256(primary_path.read_bytes()).hexdigest()
    primary=pd.read_csv(primary_path,low_memory=False)
    folder=OUT/'primary_vs_plain_conditions'
    folder.mkdir(exist_ok=True)
    results={'primary_vs_plain_conditions':comparisons(primary,
        [(name,'conditions') for name in ['composition','surface','esm650']],folder,family=True)}
    full,_=prepare_cohort(load_measurements(raw_path))
    clusters=pd.read_csv(ROOT/'reports/family_benchmark/clusters.csv',usecols=['sequence_id','identity30'])
    full=full.merge(clusters,on='sequence_id',validate='many_to_one')
    surfaces=pd.read_csv(ROOT/'reports/surface/features.csv')
    surfaces=surfaces[['sequence_id']+[c for c in surfaces if c.startswith('surface_')]]
    for variant,seed in [('augmented_metadata',None),('augmented_metadata_split47',47)]:
        assigned=sequence_publication_folds(full,random_state=seed)
        rows=align_cohort(annotate_metadata(add_amino_acid_composition(assigned)),surfaces)
        # Specifications may refer to an unused ESM set, but no embeddings enter these fits.
        results[variant]=familiar_task(rows,assigned,variant=variant,metadata=META_FEATURES,
            names=['family_conditions','composition','surface'])
        (OUT/'robustness_results.json').write_text(json.dumps(results,indent=2))
    results['checks']={'primary_predictions_unchanged':hashlib.sha256(primary_path.read_bytes()).hexdigest()==primary_hash,
        'raw_sha256_after':hashlib.sha256(raw_path.read_bytes()).hexdigest()}
    assert results['checks']['primary_predictions_unchanged']
    assert results['checks']['raw_sha256_after']==expected
    (OUT/'robustness_results.json').write_text(json.dumps(results,indent=2))
    print('ROBUSTNESS COMPLETE',json.dumps(results),flush=True)


def run_repeated_stability():
    from .metadata_diagnostics import annotate_metadata, META_FEATURES
    assert (OUT/'REPEATED_SPLIT_AMENDMENT.md').exists()
    folder=OUT/'repeated_split_stability'
    folder.mkdir(exist_ok=True)
    expected=json.loads((OUT/'analysis_plan.json').read_text())['raw_sha256']
    raw_path=ROOT/'data/raw/measurements.csv'
    assert hashlib.sha256(raw_path.read_bytes()).hexdigest()==expected
    oldpaths=[OUT/v/'predictions.csv' for v in ['familiar_family','augmented_metadata','augmented_metadata_split47','known_enzyme']]
    oldhashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in oldpaths}
    predictions={name:pd.read_csv(OUT/name/'predictions.csv',low_memory=False)
        for name in ['augmented_metadata','augmented_metadata_split47']}
    shared=set(predictions['augmented_metadata'].measurement_id)&set(predictions['augmented_metadata_split47'].measurement_id)
    results={'shared_measurement_ids':len(shared),'shared_rows_scores':{},'assignments':{}}
    for name,predicted in predictions.items():
        subfolder=folder/(name+'_shared_rows')
        subfolder.mkdir(exist_ok=True)
        common=predicted.loc[predicted.measurement_id.isin(shared)].copy()
        common['cluster']=common.identity30.map(uncertainty_groups(common,'identity30'))
        results['shared_rows_scores'][name]=comparisons(common,[('composition','family_conditions'),
            ('surface','family_conditions')],subfolder,family=True)
    full,_=prepare_cohort(load_measurements(raw_path))
    clusters=pd.read_csv(ROOT/'reports/family_benchmark/clusters.csv',usecols=['sequence_id','identity30'])
    full=full.merge(clusters,on='sequence_id',validate='many_to_one')
    surfaces=pd.read_csv(ROOT/'reports/surface/features.csv')
    surfaces=surfaces[['sequence_id']+[c for c in surfaces if c.startswith('surface_')]]
    for name,seed in [('augmented_metadata',None),('augmented_metadata_split47',47),
                      ('augmented_metadata_split17',17),('augmented_metadata_split23',23),
                      ('augmented_metadata_split71',71),('augmented_metadata_split101',101)]:
        assigned=sequence_publication_folds(full,random_state=seed)
        rows=align_cohort(annotate_metadata(add_amino_acid_composition(assigned)),surfaces)
        if name in predictions:
            baselinevariant=name+'_plain_baseline'
            familiar_task(rows,assigned,variant=baselinevariant,metadata=META_FEATURES,names=['conditions'])
            baseline=pd.read_csv(OUT/baselinevariant/'predictions.csv',low_memory=False)
            predicted=predictions[name]
            assert predicted.measurement_id.tolist()==baseline.measurement_id.tolist()
            predicted=predicted.copy()
            predicted['conditions']=baseline.conditions.to_numpy()
        else:
            familiar_task(rows,assigned,variant=name,metadata=META_FEATURES,
                names=['conditions','family_conditions','composition'])
            predicted=pd.read_csv(OUT/name/'predictions.csv',low_memory=False)
        predictions[name]=predicted
        subfolder=folder/name
        subfolder.mkdir(exist_ok=True)
        predicted.to_csv(subfolder/'predictions.csv',index=False)
        results['assignments'][name]=comparisons(predicted,[('composition','family_conditions'),
            ('composition','conditions'),('family_conditions','conditions')],subfolder,family=True)
        (folder/'results.json').write_text(json.dumps(results,indent=2))
    # Aggregate scores, never repeated rows or CV runs as independent test units.
    allrows=pd.concat(predictions.values(),ignore_index=True)
    familygroups=uncertainty_groups(allrows,'identity30')
    seqfamily=allrows[['sequence_id','identity30']].drop_duplicates()
    results['aggregate']={}
    for reference in ['family_conditions','conditions']:
        label='composition_vs_'+reference
        tables=[pd.read_csv(folder/name/(label+'_scores.csv')) for name in predictions]
        combined=aggregate_repeated_scores(tables).merge(seqfamily,on='sequence_id',validate='one_to_one')
        combined['cluster']=combined.identity30.map(familygroups)
        combined['delta']=combined.candidate_score-combined.reference_score
        stats=paired_cluster_ci(combined.rename(columns={'candidate_score':'composition','reference_score':'baseline'}))
        stats.update(candidate_mean=float(combined.candidate_score.mean()),reference_mean=float(combined.reference_score.mean()),
            min_available_runs=int(combined.available_runs.min()),max_available_runs=int(combined.available_runs.max()))
        results['aggregate'][label]=stats
        combined.to_csv(folder/(label+'_aggregate_scores.csv'),index=False)
        combined.groupby('identity30').agg(sequences=('sequence_id','size'),delta=('delta','mean'),
            candidate=('candidate_score','mean'),reference=('reference_score','mean')).to_csv(folder/(label+'_aggregate_per_family.csv'))
        pd.DataFrame([dict(omitted_family=f,delta=float(combined.loc[combined.identity30.ne(f),'delta'].mean()))
            for f in sorted(combined.identity30.unique())]).to_csv(folder/(label+'_aggregate_leave_one_family_out_scores_only.csv'),index=False)
        print('AGGREGATE',label,json.dumps(stats),flush=True)
    results['checks']=dict(raw_sha256_after=hashlib.sha256(raw_path.read_bytes()).hexdigest(),
        previous_predictions_unchanged=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in oldhashes.items()))
    assert results['checks']['raw_sha256_after']==expected
    assert results['checks']['previous_predictions_unchanged']
    (folder/'results.json').write_text(json.dumps(results,indent=2))
    print('REPEATED STABILITY COMPLETE',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['primary','robustness','repeated'],default='primary')
    options=parser.parse_args()
    {'primary':run,'robustness':run_robustness,'repeated':run_repeated_stability}[options.phase]()
