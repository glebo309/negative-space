"""Supplementary fixed-regularization learner check; no held-out parameter tuning."""
import hashlib
import inspect
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error,r2_score
from threadpoolctl import threadpool_limits

from .benchmark import prepare_cohort
from .data import load_measurements
from .embeddings import join_embeddings
from .features import add_amino_acid_composition
from .interaction_model import InteractionFeatures
from .run_benchmark import CONDITIONS
from .surface_benchmark import align_cohort,feature_sets,evaluate

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/surface_interactions'


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    raw=ROOT/'data/raw/measurements.csv'
    reference=json.loads((ROOT/'reports/family_benchmark/results.json').read_text())
    raw_hash=hashlib.sha256(raw.read_bytes()).hexdigest()
    assert raw_hash==reference['audit']['raw_sha256']
    frame,_=prepare_cohort(load_measurements(raw))
    frame=add_amino_acid_composition(frame)
    embedded=np.load(ROOT/'reports/esm650m/embeddings.npz',allow_pickle=False)
    sequences=frame[['sequence_id','sequence']].drop_duplicates().sort_values('sequence_id')
    assert np.array_equal(embedded['ids'],sequences.sequence_id)
    assert np.array_equal(embedded['hashes'],[hashlib.sha256(s.encode()).hexdigest() for s in sequences.sequence])
    frame=join_embeddings(frame,embedded['ids'].tolist(),embedded['vectors'])
    surfaces=pd.read_csv(ROOT/'reports/surface/features.csv')
    surfaces=surfaces[['sequence_id']+[c for c in surfaces if c.startswith('surface_')]]
    plan=dict(raw_sha256=raw_hash,alpha=100,solver='lsqr',tolerance=1e-6,max_iter=2000,
        esm_pca_components=16,protein_scaling='Unique training proteins only; ESM PCA also train-only',
        context='Training-only imputation/scaling and solvent one-hot encoding',
        basis='Conditions + solvent + solvent*conditions; candidates add protein + protein*solvent + protein*conditions',
        target='log1p(percent), predictions floored at log value zero before conversion',
        scope='Supplementary learner check chosen after primary surface-tree results; fixed alpha, no held-out tuning')
    plan_path=OUT/'analysis_plan.json'
    if plan_path.exists():
        assert json.loads(plan_path.read_text())==plan
    else:
        plan_path.write_text(json.dumps(plan,indent=2))
    result_path=OUT/'results.json'
    results=json.loads(result_path.read_text()) if result_path.exists() else {}
    for split in ['identity30_publication','identity30','identity40','sequence']:
        folder=OUT/split
        folder.mkdir(exist_ok=True)
        stored=pd.read_csv(ROOT/'reports/family_benchmark'/f'{split}_predictions.csv',low_memory=False)
        assert frame.measurement_id.tolist()==stored.measurement_id.tolist()
        rows=frame.copy()
        rows['fold']=stored.fold.to_numpy()
        rows['cluster']=stored.cluster.to_numpy()
        rows=align_cohort(rows,surfaces)
        assert rows.groupby('cluster').fold.nunique().max()==1
        assert rows.groupby('sequence_id').fold.nunique().max()==1
        if split.endswith('publication'):
            assert rows.groupby('doi').fold.nunique().max()==1
        sets=feature_sets(rows)
        evaluated=rows[[c for c in rows if not c.startswith(('aa_','esm_','surface_')) and c!='sequence']].copy()
        details=dict(rows=len(rows),sequences=int(rows.sequence_id.nunique()),models={},comparisons={})
        for name in ['conditions','composition','surface','esm650']:
            protein=sets[name]
            keys=['measurement_id','sequence_id','fold','cluster','doi','solvent_name','target_log1p']+CONDITIONS+protein
            signature=hashlib.sha256(pd.util.hash_pandas_object(rows[keys],index=False).values.tobytes()
                +json.dumps(plan,sort_keys=True).encode()+inspect.getsource(InteractionFeatures).encode()
                +inspect.getsource(run).encode()).hexdigest()
            cache=folder/f'{name}.npz'
            prediction=np.full(len(rows),np.nan)
            iterations=[]
            clipped=0
            if cache.exists():
                saved=np.load(cache,allow_pickle=False)
                assert saved['signature'].item()==signature
                prediction=saved['predictions']
                iterations=saved['iterations'].tolist()
                clipped=int(saved['clipped'])
            else:
                for fold in sorted(rows.fold.unique()):
                    test=rows.fold.eq(fold).to_numpy()
                    with threadpool_limits(limits=4):
                        encoder=InteractionFeatures(CONDITIONS,protein,pca_components=16 if name=='esm650' else None)
                        train_x=encoder.fit_transform(rows.loc[~test])
                        test_x=encoder.transform(rows.loc[test])
                        model=Ridge(alpha=100,solver='lsqr',tol=1e-6,max_iter=2000)
                        model.fit(train_x,rows.loc[~test,'target_log1p'])
                        raw_prediction=model.predict(test_x)
                    clipped+=int((raw_prediction<0).sum())
                    prediction[test]=np.expm1(np.maximum(raw_prediction,0))
                    iterations.append(int(np.max(model.n_iter_)))
                    print('INTERACTION',split,name,'fold',int(fold)+1,'iterations',iterations[-1],flush=True)
                assert np.isfinite(prediction).all()
                assert max(iterations)<2000
                np.savez_compressed(cache,signature=np.array(signature),predictions=prediction,
                                    iterations=np.array(iterations),clipped=np.array(clipped))
            evaluated[name]=prediction
            details['models'][name]=dict(mae=float(mean_absolute_error(rows.measured_value,prediction)),
                r2_percent=float(r2_score(rows.measured_value,prediction)),
                r2_log1p=float(r2_score(rows.target_log1p,np.log1p(prediction))),
                iterations=iterations,clipped_negative_log_predictions=clipped)
            evaluated.to_csv(folder/'predictions.csv',index=False)
        for candidate,baseline in [('composition','conditions'),('surface','conditions'),
                                    ('esm650','conditions'),('surface','composition')]:
            label=candidate+'_vs_'+baseline
            details['comparisons'][label]={}
            for kind,matched in [('matched_solvents',True),('all_conditions',False)]:
                stats,scores=evaluate(evaluated,candidate,baseline,matched=matched)
                details['comparisons'][label][kind]=stats
                scores.to_csv(folder/f'{label}_{kind}_scores.csv',index=False)
            print('INTERACTION RESULT',split,label,
                  json.dumps(details['comparisons'][label]['matched_solvents']),flush=True)
        results[split]=details
        result_path.write_text(json.dumps(results,indent=2))


if __name__=='__main__':
    run()
