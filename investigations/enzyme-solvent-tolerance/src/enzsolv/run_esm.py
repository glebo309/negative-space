"""Extract frozen ESM2 features and evaluate against saved baseline predictions."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmModel
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from .benchmark import prepare_cohort, rank_scores, paired_cluster_ci
from .data import load_measurements
from .embeddings import residue_mean, sequence_chunks, join_embeddings
from .run_benchmark import CONDITIONS

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'reports/family_benchmark'


def model_spec(size):
    models = {
        '35m': ('facebook/esm2_t12_35M_UR50D','6fbf070e65b0b7291e7bbcd451118c216cff79d8'),
        '650m': ('facebook/esm2_t33_650M_UR50D','08e4846e537177426273712802403f7ba8261b6c'),
    }
    model_id, revision = models[size]
    return ROOT/f'reports/esm{size}', model_id, revision


def extract(sequences, size='35m'):
    out, model_id, revision = model_spec(size)
    cache = out/'embeddings.npz'
    hashes = np.array([hashlib.sha256(s.encode()).hexdigest() for s in sequences.sequence])
    ids = sequences.sequence_id.to_numpy(dtype=str)
    if cache.exists():
        saved = np.load(cache,allow_pickle=False)
        assert np.array_equal(saved['ids'],ids)
        assert np.array_equal(saved['hashes'],hashes)
        assert saved['revision'].item()==revision
        return saved['vectors'], json.loads((out/'extraction.json').read_text())
    device = 'mps' if torch.backends.mps.is_available() else 'cpu'
    tokenizer = AutoTokenizer.from_pretrained(model_id,revision=revision,local_files_only=True)
    model = EsmModel.from_pretrained(model_id,revision=revision,local_files_only=True,
        use_safetensors=True,add_pooling_layer=False).eval()
    model.requires_grad_(False)
    torch.set_num_threads(4)

    def embed(chunk, destination):
        batch = tokenizer(chunk,return_tensors='pt',return_special_tokens_mask=True,truncation=False)
        special = batch.pop('special_tokens_mask').to(destination)
        batch = {k:v.to(destination) for k,v in batch.items()}
        with torch.inference_mode():
            hidden = model(**batch).last_hidden_state
            return residue_mean(hidden,batch['attention_mask'],special).cpu().numpy()[0]

    # Verify the real model agrees on CPU and GPU before processing the cohort.
    chunk = sequences.sequence.iloc[0][:128]
    cpu = embed(chunk,'cpu')
    model.to(device)
    gpu = embed(chunk,device)
    difference = float(np.max(np.abs(cpu-gpu)))
    assert np.allclose(cpu,gpu,atol=2e-3,rtol=2e-3), difference
    print('ESM loaded',device,'CPU/GPU max difference',difference,flush=True)
    start = time.monotonic()
    vectors = []
    long_sequences = []
    for i,row in enumerate(sequences.itertuples(),1):
        chunks = sequence_chunks(row.sequence)
        if len(chunks)>1:
            long_sequences.append(row.sequence_id)
        embeddings = [embed(c,device) for c in chunks]
        vectors.append(np.average(embeddings,axis=0,weights=[len(c) for c in chunks]))
        if i%25==0 or i==len(sequences):
            print('Embedded',i,'/',len(sequences),'elapsed',round(time.monotonic()-start,1),'s',flush=True)
    vectors = np.asarray(vectors,dtype=np.float32)
    assert np.isfinite(vectors).all()
    np.savez_compressed(cache,ids=ids,hashes=hashes,vectors=vectors,revision=np.array(revision))
    metadata = dict(model=model_id,revision=revision,device=device,shape=list(vectors.shape),
        seconds=time.monotonic()-start,cpu_gpu_max_difference=difference,
        long_sequence_ids=long_sequences,chunk_size=1022,
        pooling='Final layer mean over residues, excluding special and padding tokens; length-weighted nonoverlapping chunks for long proteins',
        versions={n:importlib.metadata.version(n) for n in ['torch','transformers','safetensors']})
    (out/'extraction.json').write_text(json.dumps(metadata,indent=2))
    return vectors,metadata


def run(size='35m'):
    out, _, _ = model_spec(size)
    out.mkdir(parents=True,exist_ok=True)
    frame,audit = prepare_cohort(load_measurements(ROOT/'data/raw/measurements.csv'))
    reference = json.loads((BASE/'results.json').read_text())
    raw_hash = hashlib.sha256((ROOT/'data/raw/measurements.csv').read_bytes()).hexdigest()
    assert raw_hash==reference['audit']['raw_sha256']
    sequences = frame[['sequence_id','sequence']].drop_duplicates().sort_values('sequence_id')
    # The legacy cluster export used 'sequence' as both raw sequence and split name.
    # Recover the amino-acid strings from the immutable raw cohort, verifying FASTA.
    fasta = (BASE/'sequences.fasta').read_text().strip().splitlines()
    fasta_map = dict(zip([s[1:] for s in fasta[::2]],fasta[1::2]))
    assert all(fasta_map[r.sequence_id]==r.sequence for r in sequences.itertuples())
    sequences.to_csv(out/'source_sequences.csv',index=False)
    vectors,metadata = extract(sequences,size=size)
    summary = dict(extraction=metadata,audit=audit,raw_sha256=raw_hash,
        model_settings=dict(n_estimators=250,min_samples_leaf=5,max_features=1.,random_state=11,target='log1p(percent)'),results={})
    featured = join_embeddings(frame,sequences.sequence_id.tolist(),vectors)
    columns = CONDITIONS+[c for c in featured.columns if c.startswith('esm_')]
    for split in ['identity30','identity30_publication','identity40','sequence']:
        stored = pd.read_csv(BASE/f'{split}_predictions.csv',low_memory=False)
        assert len(stored)==len(featured)
        for c in ['measurement_id','sequence_id','doi']:
            assert stored[c].astype(str).tolist()==featured[c].astype(str).tolist()
        assert np.allclose(stored.measured_value,featured.measured_value)
        for c in CONDITIONS:
            assert np.allclose(stored[c].astype(float),featured[c].astype(float),equal_nan=True)
        predicted = stored.copy()
        for fold in sorted(stored.fold.unique()):
            train = stored.fold.ne(fold).to_numpy()
            test = ~train
            assert set(stored.loc[train,'cluster']).isdisjoint(stored.loc[test,'cluster'])
            if split.endswith('publication'):
                assert set(stored.loc[train,'doi']).isdisjoint(stored.loc[test,'doi'])
            model = Pipeline([
                ('prep',ColumnTransformer([
                    ('num',SimpleImputer(strategy='median',add_indicator=True,keep_empty_features=True),columns),
                    ('cat',OneHotEncoder(handle_unknown='ignore',sparse_output=False),['solvent_name'])])),
                ('model',ExtraTreesRegressor(n_estimators=250,min_samples_leaf=5,
                    max_features=1.,n_jobs=4,random_state=11))])
            model.fit(featured.loc[train],featured.loc[train,'target_log1p'])
            predicted.loc[test,'esm'] = np.expm1(model.predict(featured.loc[test]))
            print(split,'fold',fold+1,'complete',flush=True)
        predicted.to_csv(out/f'{split}_predictions.csv',index=False)
        metrics = dict(mae_percentage_points=float(mean_absolute_error(predicted.measured_value,predicted.esm)),
            r2_percent=float(r2_score(predicted.measured_value,predicted.esm)),
            r2_log1p=float(r2_score(np.log1p(predicted.measured_value),np.log1p(predicted.esm))))
        for kind,matched in [('all_conditions',False),('matched_solvents',True)]:
            comparison = predicted.copy()
            # Existing scoring API calls the candidate column composition.
            comparison['composition'] = comparison.esm
            scores = rank_scores(comparison,matched=matched)
            statistics = paired_cluster_ci(scores)
            statistics.update(baseline=float(scores.baseline.mean()),esm=float(scores.composition.mean()),
                n_series=int(scores.n_series.sum()))
            metrics[kind] = statistics
            scores.rename(columns={'composition':'esm'}).to_csv(out/f'{split}_{kind}_scores.csv',index=False)
        summary['results'][split] = metrics
        (out/'results.json').write_text(json.dumps(summary,indent=2))
        print('RESULT',split,json.dumps(metrics),flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--size',choices=['35m','650m'],default='35m')
    run(size=parser.parse_args().size)
