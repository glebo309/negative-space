"""Remove family identity from the composition candidate on all six fixed assignments."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .benchmark import prepare_cohort,paired_cluster_ci
from .composition_ablation import assign_saved_folds,composition_model_spec
from .data import load_measurements
from .decision_utility import utility_scores
from .easier_tasks import aggregate_repeated_scores,uncertainty_groups,assert_split
from .features import add_amino_acid_composition
from .metadata_diagnostics import annotate_metadata
from .run_easier_tasks import fit_prediction
from .surface_benchmark import align_cohort,evaluate

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'reports/decision_diagnostics/easier_tasks'
OUT=ROOT/'reports/decision_diagnostics/composition_ablation'


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    raw=ROOT/'data/raw/measurements.csv'
    digest=hashlib.sha256(raw.read_bytes()).hexdigest()
    assert digest==json.loads((SOURCE/'analysis_plan.json').read_text())['raw_sha256']
    completed=json.loads((SOURCE/'repeated_split_stability/results.json').read_text())
    names=list(completed['assignments'])
    plan=dict(raw_sha256=digest,reason='Final independent reviewer identified candidate still carried family identity, a feature that harmed baseline',
        change='Composition candidate now uses exactly metadata-only baseline inputs plus 20 amino-acid frequencies; family identity omitted from both',
        fixed='Same six assignments, same supported evaluation rows, same full covered training cohort, same250trees seed11 minleaf5',
        assignments=names,scope='Bounded ablation before stop-versus-narrow decision, not heldout hyperparameter tuning',
        aggregation='Sequence-average over available assignments; count each sequence once; family/publication uncertainty components',
        positive_results='Report every assignment, primary rank metric and secondary utility, regardless of outcome')
    planpath=OUT/'analysis_plan.json'
    if planpath.exists():
        assert json.loads(planpath.read_text())==plan
    else:
        planpath.write_text(json.dumps(plan,indent=2))
    frame,_=prepare_cohort(load_measurements(raw))
    frame=annotate_metadata(add_amino_acid_composition(frame))
    surface=pd.read_csv(ROOT/'reports/surface/features.csv',usecols=['sequence_id'])
    frame=align_cohort(frame,surface)
    numeric,categorical=composition_model_spec(frame)
    results={'assignments':{}}
    allscores=[]
    allutility=[]
    allrows=[]
    for name in names:
        folder=OUT/name
        folder.mkdir(parents=True,exist_ok=True)
        assignments=pd.read_csv(SOURCE/name/'fold_assignments.csv')
        rows=assign_saved_folds(frame,assignments)
        old=pd.read_csv(SOURCE/'repeated_split_stability'/name/'predictions.csv',low_memory=False)
        predicted=old.rename(columns={'composition':'composition_with_family'}).copy()
        predicted['composition_no_family']=np.nan
        for fold in sorted(rows.fold.unique()):
            train=rows.loc[rows.fold.ne(fold)]
            ids=predicted.loc[predicted.fold.eq(fold),'measurement_id']
            test=rows.set_index('measurement_id',drop=False).loc[ids].reset_index(drop=True)
            assert_split(train,test)
            assert np.allclose(test.measured_value,predicted.loc[predicted.fold.eq(fold),'measured_value'])
            values=fit_prediction(train,test,numeric,categorical,folder/f'composition_no_family_fold{fold}.npz')
            predicted.loc[predicted.fold.eq(fold),'composition_no_family']=values
        assert np.isfinite(predicted.composition_no_family).all()
        predicted.to_csv(folder/'predictions.csv',index=False)
        stats,scores=evaluate(predicted,'composition_no_family','conditions',matched=True)
        scores.to_csv(folder/'scores.csv',index=False)
        utility=utility_scores(predicted,'composition_no_family','conditions')
        utility.to_csv(folder/'utility_scores.csv',index=False)
        allscores.append(scores)
        allutility.append(utility.assign(assignment=name))
        allrows.append(predicted)
        results['assignments'][name]=stats
        (OUT/'results.json').write_text(json.dumps(results,indent=2))
        print('NO FAMILY COMPOSITION',name,json.dumps(stats),flush=True)
    union=pd.concat(allrows,ignore_index=True)
    seqfamily=union[['sequence_id','identity30']].drop_duplicates()
    groups=uncertainty_groups(union,'identity30')
    combined=aggregate_repeated_scores(allscores).merge(seqfamily,on='sequence_id',validate='one_to_one')
    combined['cluster']=combined.identity30.map(groups)
    combined.to_csv(OUT/'aggregate_scores.csv',index=False)
    stats=paired_cluster_ci(combined.rename(columns={'candidate_score':'composition','reference_score':'baseline'}))
    stats.update(candidate_mean=float(combined.candidate_score.mean()),reference_mean=float(combined.reference_score.mean()))
    results['aggregate']=stats
    combined['delta']=combined.candidate_score-combined.reference_score
    combined.groupby('identity30').agg(sequences=('sequence_id','size'),delta=('delta','mean')).to_csv(OUT/'per_family.csv')
    utility=pd.concat(allutility,ignore_index=True)
    metrics=[c for c in utility if c.startswith(('candidate_','reference_'))]
    aggregate=utility.groupby('sequence_id',as_index=False)[metrics].mean().merge(seqfamily,on='sequence_id',validate='one_to_one')
    aggregate['cluster']=aggregate.identity30.map(groups)
    aggregate.to_csv(OUT/'aggregate_utility_scores.csv',index=False)
    results['utility']={}
    for metric in ['pairwise','top_rank','best_hit']:
        selected=aggregate[['sequence_id','cluster','candidate_'+metric,'reference_'+metric]].rename(
            columns={'candidate_'+metric:'composition','reference_'+metric:'baseline'})
        detail=paired_cluster_ci(selected)
        detail.update(candidate_mean=float(selected.composition.mean()),reference_mean=float(selected.baseline.mean()))
        results['utility'][metric]=detail
    results['raw_sha256_after']=hashlib.sha256(raw.read_bytes()).hexdigest()
    assert results['raw_sha256_after']==digest
    (OUT/'results.json').write_text(json.dumps(results,indent=2))
    print('NO FAMILY AGGREGATE',json.dumps(results),flush=True)


if __name__=='__main__':
    run()
