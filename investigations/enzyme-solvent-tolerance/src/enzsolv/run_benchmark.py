"""Run locally: PYTHONPATH=src python3 -m enzsolv.run_benchmark"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import importlib.metadata

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from .data import load_measurements
from .features import AMINO_ACIDS, add_amino_acid_composition
from .benchmark import prepare_cohort, connected_groups, paired_cluster_ci, rank_scores, assignment_table

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'reports'/'family_benchmark'
CONDITIONS = ['solvent_volume','incubation_temperature_c','assay_temperature_c',
              'incubation_ph','assay_ph','incubation_hours','solvent_logp']


def resolve_mmseqs_binary(root=ROOT):
    override = os.environ.get('MMSEQS_BIN')
    if override:
        binary = Path(override).expanduser()
        if binary.is_file():
            return binary
        raise FileNotFoundError(f'MMSEQS_BIN does not point to a file: {binary}')

    bundled = root/'tools/mmseqs/bin/mmseqs'
    if bundled.is_file():
        return bundled

    installed = shutil.which('mmseqs')
    if installed:
        return Path(installed)

    raise FileNotFoundError(
        'MMseqs2 was not found. Install mmseqs on PATH or set MMSEQS_BIN to its executable.'
    )


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    raw_path = ROOT/'data/raw/measurements.csv'
    frame, audit = prepare_cohort(load_measurements(raw_path))
    frame = add_amino_acid_composition(frame)
    audit['numeric_condition_coverage'] = {c:float(frame[c].notna().mean()) for c in CONDITIONS}
    audit['raw_sha256'] = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    audit['python'] = sys.version
    audit['versions'] = {name:importlib.metadata.version(name) for name in ['pandas','numpy','scipy','scikit-learn']}
    sequences = frame[['sequence_id','sequence']].drop_duplicates().sort_values('sequence_id')
    fasta = OUT/'sequences.fasta'
    fasta.write_text(''.join(f'>{row.sequence_id}\n{row.sequence}\n' for row in sequences.itertuples()))
    binary = resolve_mmseqs_binary(ROOT)
    audit['mmseqs_version'] = subprocess.check_output([str(binary),'version'],text=True).strip()
    hits_path = OUT/'all_pairs.tsv'
    command = [str(binary),'easy-search',str(fasta),str(fasta),str(hits_path),str(OUT/'mmseqs_tmp'),
        '--prefilter-mode','2','--alignment-mode','3','--seq-id-mode','0',
        '--min-seq-id','0.3','-c','0.8','--cov-mode','0','-e','0.001',
        '--max-seqs',str(len(sequences)),'--threads','4','-v','1',
        '--format-output','query,target,fident,qcov,tcov,evalue']
    audit['mmseqs_command'] = command
    with (OUT/'mmseqs.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    hits = pd.read_csv(hits_path,sep='\t',names=['query','target','identity','qcov','tcov','evalue'])
    ids = sequences.sequence_id.tolist()
    summary = {}
    groupings = {}
    for split in ['sequence','identity40','identity30','identity30_publication']:
        edges = []
        if split!='sequence':
            threshold = .4 if split=='identity40' else .3
            valid = hits[(hits.identity>=threshold)&(hits.qcov>=.8)&(hits.tcov>=.8)&(hits.evalue<=.001)]
            edges = list(valid[['query','target']].itertuples(index=False,name=None))
        if split.endswith('publication'):
            for _, group in frame.groupby('doi'):
                members = group.sequence_id.unique()
                edges += [(members[0],member) for member in members[1:]]
        assignments = connected_groups(ids,edges)
        groupings[split] = assignments
        evaluated = frame.copy()
        evaluated['cluster'] = evaluated.sequence_id.map(assignments)
        counts = evaluated.groupby('cluster').sequence_id.nunique()
        details = dict(n_clusters=len(counts),largest_cluster_sequences=int(counts.max()),
                       similarity_edges=len(edges))
        evaluated['fold'] = -1
        folds = list(GroupKFold(n_splits=5).split(evaluated,groups=evaluated.cluster))
        for fold,(train,test) in enumerate(folds):
            assert set(evaluated.iloc[train].cluster).isdisjoint(evaluated.iloc[test].cluster)
            if split.endswith('publication'):
                assert set(evaluated.iloc[train].doi).isdisjoint(evaluated.iloc[test].doi)
            evaluated.loc[test,'fold'] = fold
        fold_by_sequence = evaluated.groupby('sequence_id').fold.first().to_dict()
        assert all(fold_by_sequence[a]==fold_by_sequence[b] for a,b in edges)
        details['cross_fold_similarity_edges'] = 0
        details['models'] = {}
        for name in ['baseline','composition']:
            numeric = CONDITIONS + ([f'aa_{aa}' for aa in AMINO_ACIDS] if name=='composition' else [])
            model = Pipeline([
                ('prep',ColumnTransformer([
                    ('num',SimpleImputer(strategy='median',add_indicator=True,keep_empty_features=True),numeric),
                    ('cat',OneHotEncoder(handle_unknown='ignore',sparse_output=False),['solvent_name'])])),
                ('model',ExtraTreesRegressor(n_estimators=250,min_samples_leaf=5,
                    max_features=1.,n_jobs=4,random_state=11))])
            prediction = np.empty(len(evaluated))
            for train,test in folds:
                model.fit(evaluated.iloc[train],evaluated.iloc[train].target_log1p)
                prediction[test] = model.predict(evaluated.iloc[test])
            evaluated[name] = np.expm1(prediction)
            details['models'][name] = {
                'mae_percentage_points':float(mean_absolute_error(evaluated.measured_value,evaluated[name])),
                'r2_percent':float(r2_score(evaluated.measured_value,evaluated[name])),
                'r2_log1p':float(r2_score(evaluated.target_log1p,prediction))}
        for kind,matched in [('all_conditions',False),('matched_solvents',True)]:
            scores = rank_scores(evaluated,matched=matched)
            scores.to_csv(OUT/f'{split}_{kind}_scores.csv',index=False)
            details[kind] = paired_cluster_ci(scores)
            details[kind].update({name:float(scores[name].mean()) for name in ['baseline','composition']})
            details[kind]['n_series'] = int(scores.n_series.sum())
        evaluated.drop(columns=['sequence']).to_csv(OUT/f'{split}_predictions.csv',index=False)
        summary[split] = details
        print(split,json.dumps(details),flush=True)
    assignment_table(sequences,groupings).to_csv(OUT/'clusters.csv',index=False)
    payload = {'audit':audit,'settings':dict(n_estimators=250,min_samples_leaf=5,max_features=1,
        seed=11,folds=5,target='log1p(measured percent)',bootstrap='paired cluster resampling, 5000 draws'),
        'results':summary}
    (OUT/'results.json').write_text(json.dumps(payload,indent=2))
    print('AUDIT',json.dumps(audit),flush=True)


if __name__=='__main__':
    run()
