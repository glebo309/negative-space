"""Training fit and intentionally leaky row CV, never evidence of transfer."""
import hashlib
import inspect
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

from .benchmark import prepare_cohort,rank_scores
from .data import load_measurements
from .decision_diagnostics import optimistic_row_folds
from .embeddings import join_embeddings
from .run_benchmark import CONDITIONS
from .surface_benchmark import align_cohort,feature_sets,model_pipeline

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/decision_diagnostics/optimistic'


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    raw=ROOT/'data/raw/measurements.csv'
    digest=hashlib.sha256(raw.read_bytes()).hexdigest()
    assert digest==json.loads((ROOT/'reports/family_benchmark/results.json').read_text())['audit']['raw_sha256']
    frame,_=prepare_cohort(load_measurements(raw))
    saved=np.load(ROOT/'reports/esm650m/embeddings.npz',allow_pickle=False)
    sequences=frame[['sequence_id','sequence']].drop_duplicates().sort_values('sequence_id')
    assert np.array_equal(saved['ids'],sequences.sequence_id.to_numpy())
    assert np.array_equal(saved['hashes'],[hashlib.sha256(s.encode()).hexdigest() for s in sequences.sequence])
    frame=join_embeddings(frame,saved['ids'].tolist(),saved['vectors'])
    surface=pd.read_csv(ROOT/'reports/surface/features.csv')
    surface=surface[['sequence_id']+[c for c in surface if c.startswith('surface_')]]
    rows=align_cohort(frame,surface)
    sets=feature_sets(rows)
    rows['fold']=optimistic_row_folds(rows)
    rows['cluster']=rows.sequence_id
    plan=dict(raw_sha256=digest,purpose='Distinguish training/design interpolation from held-out enzyme prediction',
        modes=['five-fold shuffled row split seed11','resubstitution on complete training dataset'],
        comparators=['conditions','surface','esm650'],settings=dict(trees=250,min_samples_leaf=5,seed=11),
        inference='Descriptive diagnostics only. Rows deliberately share proteins and publications across folds; no valid transfer claim or inferential confidence interval.')
    (OUT/'analysis_plan.json').write_text(json.dumps(plan,indent=2))
    coverage=[]
    for fold in sorted(rows.fold.unique()):
        train=rows.loc[rows.fold.ne(fold)]
        test=rows.loc[rows.fold.eq(fold)]
        coverage.append(dict(fold=int(fold),rows=len(test),
            test_rows_with_known_sequence=float(test.sequence_id.isin(train.sequence_id).mean()),
            test_rows_with_known_publication=float(test.doi.isin(train.doi).mean())))
    (OUT/'overlap.json').write_text(json.dumps(coverage,indent=2))
    results={}
    for mode in ['random_rows','training_fit']:
        evaluated=rows[[c for c in rows if not c.startswith(('esm_','surface_')) and c!='sequence']].copy()
        details={}
        for name in ['conditions','surface','esm650']:
            numeric=CONDITIONS+sets[name]
            signature=hashlib.sha256(pd.util.hash_pandas_object(rows[['measurement_id','fold','target_log1p','solvent_name']+numeric],index=False).values.tobytes()+inspect.getsource(model_pipeline).encode()+mode.encode()).hexdigest()
            path=OUT/f'{mode}_{name}.npz'
            predicted=np.full(len(rows),np.nan)
            completed=[]
            if path.exists():
                cache=np.load(path,allow_pickle=False)
                assert cache['signature'].item()==signature
                predicted=cache['predictions'].copy()
                completed=cache['completed'].tolist()
            for fold in (list(sorted(rows.fold.unique())) if mode=='random_rows' else [-1]):
                if fold in completed:
                    continue
                test=rows.fold.eq(fold).to_numpy() if mode=='random_rows' else np.ones(len(rows),dtype=bool)
                train=~test if mode=='random_rows' else test
                model=model_pipeline(numeric)
                model.fit(rows.loc[train],rows.loc[train,'target_log1p'])
                predicted[test]=np.expm1(model.predict(rows.loc[test]))
                completed.append(int(fold))
                temporary=path.with_suffix('.tmp.npz')
                np.savez_compressed(temporary,signature=np.array(signature),predictions=predicted,completed=np.array(completed))
                temporary.replace(path)
                print('OPTIMISTIC',mode,name,'fold',fold,'complete',flush=True)
            assert np.isfinite(predicted).all()
            evaluated[name]=predicted
            details[name]=dict(r2_log1p=float(r2_score(rows.target_log1p,np.log1p(predicted))))
        evaluated.to_csv(OUT/f'{mode}_predictions.csv',index=False)
        for candidate in ['surface','esm650']:
            scored=evaluated.copy()
            scored[['baseline','composition']]=scored[['conditions',candidate]].to_numpy()
            scores=rank_scores(scored,matched=True)
            scores.to_csv(OUT/f'{mode}_{candidate}_scores.csv',index=False)
            details[candidate].update(mean_matched_spearman=float(scores.composition.mean()),
                delta=float((scores.composition-scores.baseline).mean()),n_sequences=len(scores),n_series=int(scores.n_series.sum()))
            details['conditions']['mean_matched_spearman']=float(scores.baseline.mean())
        results[mode]=details
        (OUT/'results.json').write_text(json.dumps(results,indent=2))
        print('OPTIMISTIC RESULT',mode,json.dumps(details),flush=True)
    assert hashlib.sha256(raw.read_bytes()).hexdigest()==digest


if __name__=='__main__':
    run()
