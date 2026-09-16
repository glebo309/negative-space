"""Direct within-experiment ranking target, predeclared exploratory sensitivity."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

from .benchmark import prepare_cohort
from .data import load_measurements
from .decision_diagnostics import matched_rank_targets,repeat_groups
from .embeddings import join_embeddings
from .features import add_amino_acid_composition
from .metadata_diagnostics import META_FEATURES,annotate_metadata,cohort_masks
from .run_benchmark import CONDITIONS
from .run_surface_benchmark import fit_oof
from .surface_benchmark import align_cohort,feature_sets,evaluate

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/decision_diagnostics'


def original_series_controls(frame):
    eligible=frame.loc[matched_rank_targets(frame).notna(),'measurement_id']
    results={}
    for split in ['identity30_publication','identity30']:
        stored=pd.read_csv(ROOT/'reports/family_benchmark'/f'{split}_predictions.csv',low_memory=False)
        for label in ['strong','moderate','no_protein_signal','shuffled']:
            path=OUT/'controls'/split/label/'predictions.csv'
            generated=pd.read_csv(path)
            # Preserve original metadata/series eligibility but use the planted target.
            rows=frame.drop(columns=['target_log1p']).merge(generated,on=['measurement_id','sequence_id'],validate='one_to_one')
            rows=rows.loc[rows.measurement_id.isin(eligible)].copy()
            rows['measured_value']=np.expm1(rows.target_log1p)
            stats,scores=evaluate(rows,'planted_protein','conditions',matched=True)
            results[split+'/'+label]=stats
    (OUT/'controls/original_eligible_series_results.json').write_text(json.dumps(results,indent=2))
    return results


def repeat_coverage(frame):
    base=['sequence_id','solvent_name','solvent_volume']
    settings={'solvent_fraction_only':base,
        'plus_duration':base+['incubation_hours'],
        'plus_incubation_temperature':base+['incubation_hours','incubation_temperature_c'],
        'plus_incubation_ph':base+['incubation_hours','incubation_temperature_c','incubation_ph']}
    result={}
    for name,keys in settings.items():
        matches=repeat_groups(frame,keys)
        matches.to_csv(OUT/'coverage'/f'{name}_relaxed_matches.csv',index=False)
        result[name]=dict(keys=keys,complete_rows=len(frame.dropna(subset=keys)),
            groups=len(matches),sequences=int(matches.sequence_id.nunique()),
            median_cross_paper_range=float(matches.range_percentage_points.median()) if len(matches) else None,
            interpretation='Partial condition matches, not certified experimental replicates; differences cannot estimate measurement noise')
    (OUT/'coverage/relaxed_summary.json').write_text(json.dumps(result,indent=2))


def run():
    raw=ROOT/'data/raw/measurements.csv'
    digest=hashlib.sha256(raw.read_bytes()).hexdigest()
    assert digest==json.loads((ROOT/'reports/family_benchmark/results.json').read_text())['audit']['raw_sha256']
    frame,audit=prepare_cohort(load_measurements(raw))
    original_series_controls(frame)
    repeat_coverage(frame)
    out=OUT/'direct_ranking'
    out.mkdir(parents=True,exist_ok=True)
    frame=annotate_metadata(add_amino_acid_composition(frame))
    frame['rank_target']=matched_rank_targets(frame)
    saved=np.load(ROOT/'reports/esm650m/embeddings.npz',allow_pickle=False)
    sequences=frame[['sequence_id','sequence']].drop_duplicates().sort_values('sequence_id')
    assert np.array_equal(saved['ids'],sequences.sequence_id.to_numpy())
    assert np.array_equal(saved['hashes'],[hashlib.sha256(s.encode()).hexdigest() for s in sequences.sequence])
    frame=join_embeddings(frame,saved['ids'].tolist(),saved['vectors'])
    surfaces=pd.read_csv(ROOT/'reports/surface/features.csv')
    surfaces=surfaces[['sequence_id']+[c for c in surfaces if c.startswith('surface_')]]
    plan=dict(raw_sha256=digest,
        purpose='Separate relative solvent ranking from heterogeneous absolute activity scales',
        target='Average replicate rows by solvent within full matched series, rank averages, normalize (midrank-1)/(n_solvents-1). Exclude <3-solvent or constant series.',
        leakage='Matched series cannot cross strict sequence/publication folds. Targets computed only from their own series; no across-series normalization.',
        models=['metadata conditions','metadata plus composition','metadata plus surface','metadata plus ESM650'],
        settings=dict(trees=250,seed=11,min_samples_leaf=5),
        cohorts=['all scoreable series','aqueous_explicit scoreable series'],
        split='original identity30_publication folds retained',
        score='Original real-activity matched-solvent Spearman. fit_oof applies a monotone expm1 transform; output is a ranking score, never activity percent.',
        status='Exploratory, specified after absolute-regression results; no heldout tuning or multiplicity correction')
    planpath=out/'analysis_plan.json'
    if planpath.exists():
        assert json.loads(planpath.read_text())==plan
    else:
        planpath.write_text(json.dumps(plan,indent=2))
    stored=pd.read_csv(ROOT/'reports/family_benchmark/identity30_publication_predictions.csv',low_memory=False)
    assert frame.measurement_id.tolist()==stored.measurement_id.tolist()
    frame['fold']=stored.fold.to_numpy()
    frame['cluster']=stored.cluster.to_numpy()
    masks=cohort_masks(frame)
    results={}
    for variant in ['all','aqueous_explicit']:
        rows=align_cohort(frame.loc[masks[variant]&frame.rank_target.notna()],surfaces)
        # fit_oof's historical column name is retained only as its generic training-target interface.
        rows['original_target_log1p']=rows.target_log1p
        rows['target_log1p']=rows.rank_target
        sets=feature_sets(rows)
        folder=out/variant
        folder.mkdir(parents=True,exist_ok=True)
        evaluated=rows[[c for c in rows if not c.startswith(('aa_','esm_','surface_')) and c!='sequence']].copy()
        details=dict(rows=len(rows),sequences=int(rows.sequence_id.nunique()),publications=int(rows.doi.nunique()),comparisons={})
        for name in ['conditions','composition','surface','esm650']:
            evaluated[name]=fit_oof(rows,CONDITIONS+META_FEATURES+sets[name],folder/f'{name}.npz',publication_holdout=True)
            evaluated.to_csv(folder/'predictions.csv',index=False)
        for candidate,reference in [('composition','conditions'),('surface','conditions'),('esm650','conditions'),('surface','composition')]:
            label=candidate+'_vs_'+reference
            stats,scores=evaluate(evaluated,candidate,reference,matched=True)
            details['comparisons'][label]=stats
            scores.to_csv(folder/f'{label}_scores.csv',index=False)
            print('DIRECT RANK',variant,label,json.dumps(stats),flush=True)
        activity=pd.read_csv(OUT/'metadata'/variant/'identity30_publication/predictions.csv',low_memory=False)
        activity=activity[['measurement_id','measured_value','metadata']].rename(
            columns={'measured_value':'source_value','metadata':'activity_metadata'})
        cross=evaluated.merge(activity,on='measurement_id',validate='one_to_one')
        assert len(cross)==len(evaluated) and np.allclose(cross.measured_value,cross.source_value)
        details['versus_activity_trained_metadata']={}
        for candidate in ['composition','surface','esm650']:
            stats,scores=evaluate(cross,candidate,'activity_metadata',matched=True)
            details['versus_activity_trained_metadata'][candidate]=stats
            scores.to_csv(folder/f'{candidate}_vs_activity_metadata_scores.csv',index=False)
            print('CROSS OBJECTIVE',variant,candidate,json.dumps(stats),flush=True)
        results[variant]=details
        (out/'results.json').write_text(json.dumps(results,indent=2))
    assert hashlib.sha256(raw.read_bytes()).hexdigest()==digest


if __name__=='__main__':
    run()
