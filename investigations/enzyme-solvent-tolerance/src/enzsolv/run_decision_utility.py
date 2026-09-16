"""Score saved predictions without fitting or choosing models from these outcomes."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
from .benchmark import paired_cluster_ci
from .decision_utility import utility_scores


def run(source,candidate,reference,output):
    source,output=Path(source),Path(output)
    rows=pd.read_csv(source,low_memory=False)
    assert rows.measurement_id.is_unique
    scores=utility_scores(rows,candidate,reference)
    output.mkdir(parents=True,exist_ok=True)
    scores.to_csv(output/'scores.csv',index=False)
    results=dict(source=str(source.resolve()),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        candidate=candidate,reference=reference,
        status='Exploratory secondary decision metrics, not new independent confirmation or model-selection targets',
        definitions=dict(pairwise='Fraction of non-tied observed solvent pairs correctly ordered; prediction ties receive half credit',
            top_rank='Observed within-panel normalized midrank of top predicted solvent, 0 worst to 1 best',
            best_hit='Fraction selecting an observed best solvent. Uniform expected score across exactly tied top predictions.'),
        weighting='Average panels within sequence, then sequences; paired bootstrap on saved uncertainty groups, conditional on fits',metrics={})
    for metric in ['pairwise','top_rank','best_hit']:
        values=scores[['sequence_id','cluster',f'candidate_{metric}',f'reference_{metric}']].rename(
            columns={f'candidate_{metric}':'composition',f'reference_{metric}':'baseline'})
        stats=paired_cluster_ci(values)
        stats.update(candidate_mean=float(values.composition.mean()),reference_mean=float(values.baseline.mean()),
                     n_series=int(scores.n_series.sum()))
        results['metrics'][metric]=stats
    (output/'results.json').write_text(json.dumps(results,indent=2))
    print(json.dumps(results,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--predictions',required=True)
    parser.add_argument('--candidate',required=True)
    parser.add_argument('--reference',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    run(args.predictions,args.candidate,args.reference,args.output)
