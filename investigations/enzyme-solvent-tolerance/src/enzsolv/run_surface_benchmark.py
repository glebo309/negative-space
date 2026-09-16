"""Run fixed-fold surface ablations, then cached ESM comparisons and sensitivities."""
import argparse
import hashlib
import importlib.metadata
import inspect
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, r2_score

from .benchmark import prepare_cohort
from .data import load_measurements
from .features import add_amino_acid_composition
from .embeddings import join_embeddings
from .run_benchmark import CONDITIONS
from .surface import summarize_surface
from .surface_benchmark import align_cohort, feature_sets, model_pipeline, evaluate

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'reports/surface_benchmark'
SPLITS = ['identity30_publication','identity30','identity40','sequence']


def fit_oof(rows, numeric, cache, seed=11, trees=250, publication_holdout=False):
    if rows.groupby('cluster').fold.nunique().max()!=1:
        raise ValueError('A cluster crosses folds')
    if rows.groupby('sequence_id').fold.nunique().max()!=1:
        raise ValueError('A sequence crosses folds')
    if publication_holdout and rows.groupby('doi').fold.nunique().max()!=1:
        raise ValueError('A publication crosses folds')
    keys = ['measurement_id','sequence_id','cluster','doi','fold','solvent_name','target_log1p']+numeric
    signature = hashlib.sha256(pd.util.hash_pandas_object(rows[keys],index=False).values.tobytes()
        +json.dumps(dict(numeric=numeric,seed=seed,trees=trees,
              sklearn=importlib.metadata.version('scikit-learn'))).encode()
        +inspect.getsource(model_pipeline).encode()+inspect.getsource(fit_oof).encode()).hexdigest()
    predicted = np.full(len(rows),np.nan)
    completed = []
    if cache.exists():
        saved = np.load(cache,allow_pickle=False)
        if saved['signature'].item()!=signature:
            raise ValueError('Cache input or model signature mismatch: '+str(cache))
        predicted = saved['predictions'].copy()
        completed = saved['completed'].tolist()
    for fold in sorted(rows.fold.unique()):
        if fold in completed:
            continue
        test = rows.fold.eq(fold).to_numpy()
        train = ~test
        if not train.any() or not test.any():
            raise ValueError('Empty training or validation fold')
        model = model_pipeline(numeric,seed=seed,trees=trees)
        model.fit(rows.loc[train],rows.loc[train,'target_log1p'])
        predicted[test] = np.expm1(model.predict(rows.loc[test]))
        completed.append(int(fold))
        cache.parent.mkdir(parents=True,exist_ok=True)
        temporary = cache.with_suffix('.tmp.npz')
        np.savez_compressed(temporary,signature=np.array(signature),predictions=predicted,
                            completed=np.array(completed))
        temporary.replace(cache)
        print(cache.stem,'fold',int(fold)+1,'complete',flush=True)
    assert np.isfinite(predicted).all()
    return predicted


def load_variants():
    manifest = json.loads((ROOT/'reports/surface/manifest.json').read_text())
    accepted = [r for r in manifest if r['status']=='ok']
    variants = {key:[] for key in ['primary','rsa10','rsa30','confidence70']}
    for entry in accepted:
        sid = entry['sequence_id']
        table = pd.read_csv(ROOT/'reports/surface/residues'/f'{sid}.csv')
        for name,cutoff in [('primary',.2),('rsa10',.1),('rsa30',.3)]:
            variants[name].append(dict(sequence_id=sid,**summarize_surface(table,cutoff)))
        confident = table.loc[table.plddt.ge(70)]
        if entry['mean_plddt']>=70 and confident.rsa.ge(.2).any():
            variants['confidence70'].append(dict(sequence_id=sid,**summarize_surface(confident,.2)))
    return {k:pd.DataFrame(v) for k,v in variants.items()}


def run(phase='all'):
    OUT.mkdir(parents=True,exist_ok=True)
    raw_path = ROOT/'data/raw/measurements.csv'
    frame,audit = prepare_cohort(load_measurements(raw_path))
    raw_hash = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    baseline = json.loads((ROOT/'reports/family_benchmark/results.json').read_text())
    assert raw_hash==baseline['audit']['raw_sha256']
    frame = add_amino_acid_composition(frame)
    saved = np.load(ROOT/'reports/esm650m/embeddings.npz',allow_pickle=False)
    sequences = frame[['sequence_id','sequence']].drop_duplicates().sort_values('sequence_id')
    assert np.array_equal(saved['ids'],sequences.sequence_id.to_numpy())
    assert np.array_equal(saved['hashes'],[hashlib.sha256(s.encode()).hexdigest() for s in sequences.sequence])
    assert saved['revision'].item()=='08e4846e537177426273712802403f7ba8261b6c'
    frame = join_embeddings(frame,saved['ids'].tolist(),saved['vectors'])
    variants = load_variants()
    plan = dict(raw_sha256=raw_hash,audit=audit,
        variants={k:dict(sequences=len(v),measurements=int(frame.sequence_id.isin(v.sequence_id).sum()))
                  for k,v in variants.items()},
        primary='27 surface descriptors at RSA >=0.20 plus conditions versus conditions-only, same rows/folds',
        controls=['whole-sequence composition','ESM650','ESM650 plus surface'],
        sensitivities=['RSA >=0.10','RSA >=0.30',
            'protein mean pLDDT >=70 and descriptors restricted to residues with pLDDT >=70',
            'seeds 23 and 47 in the primary 30%-identity-plus-publication split'],
        settings=dict(trees=250,min_samples_leaf=5,max_features=1.,seed=11,
                      target='log1p(percent)',folds='unchanged saved assignments after row filtering',
                      bootstrap='5000 paired component draws, conditional on fits'),
        inference='Exploratory sensitivity analyses, no selection of a winning cutoff; intervals are not multiplicity-adjusted',
        versions={k:importlib.metadata.version(k) for k in ['numpy','pandas','scikit-learn']})
    plan_path = OUT/'analysis_plan.json'
    if plan_path.exists():
        assert json.loads(plan_path.read_text())==plan
    else:
        plan_path.write_text(json.dumps(plan,indent=2))
    results_path = OUT/'results.json'
    results = json.loads(results_path.read_text()) if results_path.exists() else {}
    jobs = []
    if phase in ['all','surface']:
        jobs += [('primary',split,11,['conditions','composition','surface']) for split in SPLITS]
    if phase in ['all','sensitivity']:
        jobs += [(variant,split,11,['conditions','composition','surface'])
                 for variant in ['rsa10','rsa30','confidence70'] for split in SPLITS]
        jobs += [('primary','identity30_publication',seed,['conditions','composition','surface'])
                 for seed in [23,47]]
    if phase in ['all','esm']:
        jobs += [('primary',split,11,['conditions','composition','surface','esm650','esm650_surface'])
                 for split in SPLITS]
    for variant,split,seed,names in jobs:
        start = time.monotonic()
        stored = pd.read_csv(ROOT/'reports/family_benchmark'/f'{split}_predictions.csv',low_memory=False)
        assert stored.measurement_id.tolist()==frame.measurement_id.tolist()
        assert stored.sequence_id.tolist()==frame.sequence_id.tolist()
        assert np.allclose(stored.measured_value,frame.measured_value)
        rows = frame.copy()
        rows['fold'] = stored.fold.to_numpy()
        rows['cluster'] = stored.cluster.to_numpy()
        rows = align_cohort(rows,variants[variant])
        sets = feature_sets(rows)
        key = f'{variant}/{split}/seed{seed}'
        folder = OUT/key
        folder.mkdir(parents=True,exist_ok=True)
        evaluated = rows[[c for c in rows if not c.startswith(('esm_','aa_','surface_')) and c!='sequence']].copy()
        details = dict(rows=len(rows),sequences=int(rows.sequence_id.nunique()),
            clusters=int(rows.cluster.nunique()),models={},comparisons={})
        for name in names:
            numeric = CONDITIONS+sets[name]
            print('Fitting',key,name,len(numeric),'numeric features',flush=True)
            evaluated[name] = fit_oof(rows,numeric,folder/f'{name}.npz',seed=seed,
                                     publication_holdout=split.endswith('publication'))
            details['models'][name] = dict(mae=float(mean_absolute_error(rows.measured_value,evaluated[name])),
                r2_percent=float(r2_score(rows.measured_value,evaluated[name])),
                r2_log1p=float(r2_score(rows.target_log1p,np.log1p(evaluated[name]))))
            evaluated.to_csv(folder/'predictions.csv',index=False)
        pairs = [(name,'conditions') for name in names if name!='conditions']
        if 'surface' in names and 'composition' in names:
            pairs.append(('surface','composition'))
        if 'esm650_surface' in names:
            pairs += [('esm650_surface','esm650'),('esm650_surface','surface')]
        for candidate,reference in pairs:
            label = candidate+'_vs_'+reference
            details['comparisons'][label] = {}
            for kind,matched in [('matched_solvents',True),('all_conditions',False)]:
                statistics,scores = evaluate(evaluated,candidate,reference,matched=matched)
                scores.to_csv(folder/f'{label}_{kind}_scores.csv',index=False)
                details['comparisons'][label][kind] = statistics
            print('RESULT',key,label,json.dumps(details['comparisons'][label]['matched_solvents']),flush=True)
        details['seconds_this_invocation'] = time.monotonic()-start
        results[key] = details
        results_path.write_text(json.dumps(results,indent=2))
    print('Completed phase',phase,flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['all','surface','sensitivity','esm'],default='all')
    run(phase=parser.parse_args().phase)
