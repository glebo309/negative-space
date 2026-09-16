"""Source-motivated metadata sensitivity with unchanged group-held-out folds."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error

from .benchmark import prepare_cohort
from .data import load_measurements
from .embeddings import join_embeddings
from .features import add_amino_acid_composition
from .metadata_diagnostics import META_FEATURES,annotate_metadata,cohort_masks
from .run_benchmark import CONDITIONS
from .run_surface_benchmark import fit_oof
from .surface_benchmark import align_cohort,feature_sets,evaluate

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/decision_diagnostics/metadata'


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    raw=ROOT/'data/raw/measurements.csv'
    digest=hashlib.sha256(raw.read_bytes()).hexdigest()
    assert digest==json.loads((ROOT/'reports/family_benchmark/results.json').read_text())['audit']['raw_sha256']
    frame,audit=prepare_cohort(load_measurements(raw))
    frame=annotate_metadata(add_amino_acid_composition(frame))
    surfaces=pd.read_csv(ROOT/'reports/surface/features.csv')
    surfaces=surfaces[['sequence_id']+[c for c in surfaces if c.startswith('surface_')]]
    saved=np.load(ROOT/'reports/esm650m/embeddings.npz',allow_pickle=False)
    sequences=frame[['sequence_id','sequence']].drop_duplicates().sort_values('sequence_id')
    assert np.array_equal(saved['ids'],sequences.sequence_id.to_numpy())
    assert np.array_equal(saved['hashes'],[hashlib.sha256(s.encode()).hexdigest() for s in sequences.sequence])
    frame=join_embeddings(frame,saved['ids'].tolist(),saved['vectors'])
    masks=cohort_masks(frame)
    counts={name:dict(rows=int(mask.sum()),sequences=int(frame.loc[mask,'sequence_id'].nunique()),
                     publications=int(frame.loc[mask,'doi'].nunique())) for name,mask in masks.items()}
    frame[['measurement_id','sequence_id','doi']+META_FEATURES].assign(
        mode_qualified=masks['mode_qualified'],aqueous_explicit=masks['aqueous_explicit']).to_csv(OUT/'inclusion_ledger.csv',index=False)
    plan=dict(raw_sha256=digest,cohort_counts_before_surface=counts,
        purpose='Test source-audit-motivated context qualification, not optimize exclusions for performance',
        cohort_rules='See metadata_diagnostics.py and source_audit/REPORT.md. Exclude final-comment assay/control uncertainty and explicit organic references; strict requires literal aqueous-phase procedure.',
        jobs=['all/identity30_publication','aqueous_explicit/identity30_publication','aqueous_explicit/identity30'],
        comparators=['conditions','metadata','composition+metadata','surface+metadata','esm650+metadata'],
        settings=dict(trees=250,seed=11,folds='original folds retained after filtering',target='log1p(percent)'),
        interpretation='Exploratory conditional bootstrap intervals, no multiplicity adjustment; absence of a curator flag is not confirmation of correct labels',
        control_strata='Score the same out-of-fold predictions separately, no subgroup-specific training',
        note='mode_qualified exported as inclusion ledger; aqueous_explicit is the prespecified smaller training sensitivity')
    planpath=OUT/'analysis_plan.json'
    if planpath.exists():
        assert json.loads(planpath.read_text())==plan
    else:
        planpath.write_text(json.dumps(plan,indent=2))
    print('METADATA PLAN',json.dumps(plan),flush=True)
    results={}
    for variant,split in [('all','identity30_publication'),('aqueous_explicit','identity30_publication'),
                           ('aqueous_explicit','identity30')]:
        stored=pd.read_csv(ROOT/'reports/family_benchmark'/f'{split}_predictions.csv',low_memory=False)
        assert frame.measurement_id.tolist()==stored.measurement_id.tolist()
        rows=frame.copy()
        rows['fold']=stored.fold.to_numpy()
        rows['cluster']=stored.cluster.to_numpy()
        rows=align_cohort(rows.loc[masks[variant]],surfaces)
        sets=feature_sets(rows)
        folder=OUT/variant/split
        folder.mkdir(parents=True,exist_ok=True)
        evaluated=rows[[c for c in rows if not c.startswith(('esm_','aa_','surface_')) and c!='sequence']].copy()
        details=dict(rows=len(rows),sequences=int(rows.sequence_id.nunique()),publications=int(rows.doi.nunique()),
                     clusters=int(rows.cluster.nunique()),models={},comparisons={},strata={})
        for name,extra in [('conditions',[]),('metadata',META_FEATURES),
                           ('composition',META_FEATURES+sets['composition']),
                           ('surface',META_FEATURES+sets['surface']),('esm650',META_FEATURES+sets['esm650'])]:
            evaluated[name]=fit_oof(rows,CONDITIONS+extra,folder/f'{name}.npz',
                                   publication_holdout=split.endswith('publication'))
            details['models'][name]=dict(mae=float(mean_absolute_error(rows.measured_value,evaluated[name])))
            evaluated.to_csv(folder/'predictions.csv',index=False)
        for candidate,reference in [('metadata','conditions'),('composition','metadata'),('surface','metadata'),
                                     ('esm650','metadata'),('surface','composition')]:
            label=candidate+'_vs_'+reference
            stats,scores=evaluate(evaluated,candidate,reference,matched=True)
            details['comparisons'][label]=stats
            scores.to_csv(folder/f'{label}_scores.csv',index=False)
            print('METADATA RESULT',variant,split,label,json.dumps(stats),flush=True)
        for kind,column in [('water','meta_control_water'),('nonincubated','meta_control_nonincubated')]:
            subset=evaluated.loc[evaluated[column].eq(1)]
            details['strata'][kind]={}
            for candidate in ['composition','surface','esm650']:
                stats,scores=evaluate(subset,candidate,'metadata',matched=True)
                details['strata'][kind][candidate]=stats
                scores.to_csv(folder/f'{kind}_{candidate}_scores.csv',index=False)
        results[variant+'/'+split]=details
        (OUT/'results.json').write_text(json.dumps(results,indent=2))
    assert hashlib.sha256(raw.read_bytes()).hexdigest()==digest
    print('METADATA COMPLETE',flush=True)


if __name__=='__main__':
    run()
