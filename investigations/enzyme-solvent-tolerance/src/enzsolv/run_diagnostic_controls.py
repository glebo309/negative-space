"""Audit available replication and test a planted signal through real CV/scoring."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

from .benchmark import prepare_cohort
from .data import load_measurements
from .decision_diagnostics import support_tables,repeat_groups,planted_target
from .run_benchmark import CONDITIONS
from .run_surface_benchmark import fit_oof
from .surface_benchmark import align_cohort,evaluate

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/decision_diagnostics'


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    frame,audit=prepare_cohort(load_measurements(ROOT/'data/raw/measurements.csv'))
    groups=pd.read_csv(ROOT/'reports/family_benchmark/clusters.csv')
    frame=frame.merge(groups[['sequence_id','identity30','identity40']],on='sequence_id',validate='many_to_one')
    coverage=OUT/'coverage'
    coverage.mkdir(exist_ok=True)
    summary={'cohort':audit}
    for family in ['identity30','identity40']:
        seq,fam=support_tables(frame,family)
        seq.to_csv(coverage/'sequence_support.csv',index=False)
        fam.sort_values('sequences',ascending=False).to_csv(coverage/f'{family}_support.csv',index=False)
        summary[family]={'families':len(fam),'families_with_5_sequences_and_5_papers':
            int((fam.sequences.ge(5)&fam.papers.ge(5)).sum()),
            'sequences_with_multiple_papers':int(seq.papers.ge(2).sum()),
            'rows_from_repeated_sequences':int(seq.loc[seq.papers.ge(2),'rows'].sum())}
    core=['sequence_id','solvent_name','solvent_volume','incubation_hours',
          'incubation_temperature_c','incubation_ph','assay_temperature_c','assay_ph']
    for name,keys in [('numeric_conditions',core),
        ('full_recorded_context',core+['assay_solution','substrates','substrate_concentrations','procedure','comments'])]:
        repeats=repeat_groups(frame,keys)
        repeats.to_csv(coverage/f'{name}_cross_paper_repeats.csv',index=False)
        summary[name]={'cross_paper_groups':len(repeats),
                      'rows_in_complete_groups':int(repeats.rows.sum()) if len(repeats) else 0}
    (coverage/'summary.json').write_text(json.dumps(summary,indent=2))
    print('COVERAGE',json.dumps(summary),flush=True)
    surface=pd.read_csv(ROOT/'reports/surface/features.csv')[['sequence_id','surface_area_acidic_fraction']]
    mean=surface.surface_area_acidic_fraction.mean()
    sd=surface.surface_area_acidic_fraction.std()
    surface['protein_signal']=(surface.surface_area_acidic_fraction-mean)/sd
    control_dir=OUT/'controls'
    control_dir.mkdir(exist_ok=True)
    plan=dict(generator='log target = 5 + 0.3*tanh(logP/2) + strength*clip(z_acidic_ASA,-2,2)*tanh(logP/2) + Gaussian noise',
        seed=11,noise=.1,signal_mean=float(mean),signal_sd=float(sd),
        controls=['strength1.5','strength0.5','strength0.0','shuffled_strong_protein_feature'],
        scope='Planted simple global signal on real design; passing does not prove sensitivity to arbitrary biology',
        raw_data_unchanged=True,positive_target='strong-signal protein candidate rho >0.7 and improvement >0.3 in both strict splits')
    (control_dir/'analysis_plan.json').write_text(json.dumps(plan,indent=2))
    results={}
    for split in ['identity30_publication','identity30']:
        stored=pd.read_csv(ROOT/'reports/family_benchmark'/f'{split}_predictions.csv',low_memory=False)
        assert frame.measurement_id.tolist()==stored.measurement_id.tolist()
        rows=frame.copy()
        rows['fold']=stored.fold.to_numpy()
        rows['cluster']=stored.cluster.to_numpy()
        rows=align_cohort(rows,surface)
        for label,strength,shuffle in [('strong',1.5,False),('moderate',.5,False),
                                      ('no_protein_signal',0.,False),('shuffled',1.5,True)]:
            synthetic=rows.copy()
            synthetic['target_log1p']=planted_target(synthetic,strength=strength)
            synthetic['measured_value']=np.expm1(synthetic.target_log1p)
            if shuffle:
                unique=synthetic[['sequence_id','protein_signal']].drop_duplicates().set_index('sequence_id')
                permuted=np.random.default_rng(47).permutation(unique.protein_signal.to_numpy())
                mapping=dict(zip(unique.index,permuted))
                synthetic['protein_signal']=synthetic.sequence_id.map(mapping)
            folder=control_dir/split/label
            folder.mkdir(parents=True,exist_ok=True)
            for name,extra in [('conditions',[]),('planted_protein',['protein_signal'])]:
                synthetic[name]=fit_oof(synthetic,CONDITIONS+extra,folder/f'{name}.npz',
                                        publication_holdout=split.endswith('publication'))
            statistics,scores=evaluate(synthetic,'planted_protein','conditions',matched=True)
            scores.to_csv(folder/'scores.csv',index=False)
            synthetic[['measurement_id','sequence_id','fold','cluster','target_log1p',
                       'conditions','planted_protein']].to_csv(folder/'predictions.csv',index=False)
            results[split+'/'+label]=statistics
            (control_dir/'results.json').write_text(json.dumps(results,indent=2))
            print('CONTROL',split,label,json.dumps(statistics),flush=True)


if __name__=='__main__':
    run()
