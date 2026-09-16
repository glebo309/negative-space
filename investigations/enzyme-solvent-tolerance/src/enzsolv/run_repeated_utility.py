"""Aggregate secondary utility scores across the fixed six familiar-family assignments."""
import json
from pathlib import Path
import pandas as pd
from .benchmark import paired_cluster_ci
from .decision_utility import utility_scores
from .easier_tasks import uncertainty_groups

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'reports/decision_diagnostics/easier_tasks/repeated_split_stability'
OUT=ROOT/'reports/decision_diagnostics/utility/repeated_familiar'


def run():
    evidence=json.loads((SOURCE/'results.json').read_text())
    assert 'aggregate' in evidence and evidence['checks']['previous_predictions_unchanged']
    OUT.mkdir(parents=True,exist_ok=True)
    records=[]
    predictions=[]
    for name in evidence['assignments']:
        frame=pd.read_csv(SOURCE/name/'predictions.csv',low_memory=False)
        predictions.append(frame)
        scores=utility_scores(frame,'composition','conditions')
        scores['assignment']=name
        records.append(scores)
    allscores=pd.concat(records,ignore_index=True)
    allscores.to_csv(OUT/'all_assignment_scores.csv',index=False)
    allrows=pd.concat(predictions,ignore_index=True)
    seqfamily=allrows[['sequence_id','identity30']].drop_duplicates()
    assert seqfamily.sequence_id.is_unique
    groups=uncertainty_groups(allrows,'identity30')
    metrics=[c for c in allscores if c.startswith(('candidate_','reference_'))]
    combined=allscores.groupby('sequence_id',as_index=False)[metrics].mean()
    runs=allscores.groupby('sequence_id').assignment.nunique().rename('available_runs').reset_index()
    combined=combined.merge(runs,on='sequence_id',validate='one_to_one').merge(seqfamily,on='sequence_id',validate='one_to_one')
    combined['cluster']=combined.identity30.map(groups)
    combined.to_csv(OUT/'sequence_averaged_scores.csv',index=False)
    result=dict(candidate='composition+metadata+family',reference='metadata without family identity',
        assignments=list(evidence['assignments']),
        interpretation='Exploratory secondary utility. Average within sequence across available assignments; each sequence counted once. Repeated folds are not independent replication.',
        available_runs=dict(minimum=int(combined.available_runs.min()),maximum=int(combined.available_runs.max())),metrics={})
    for metric in ['pairwise','top_rank','best_hit']:
        selected=combined[['sequence_id','cluster','candidate_'+metric,'reference_'+metric]].rename(
            columns={'candidate_'+metric:'composition','reference_'+metric:'baseline'})
        stats=paired_cluster_ci(selected)
        stats.update(candidate_mean=float(selected.composition.mean()),reference_mean=float(selected.baseline.mean()))
        result['metrics'][metric]=stats
    (OUT/'results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    run()
