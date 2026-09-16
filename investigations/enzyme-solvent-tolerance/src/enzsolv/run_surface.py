"""Retrieve database structures and export sequence-verified residue accessibility."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import importlib.metadata
import json
from pathlib import Path
import time
from urllib.parse import quote

import pandas as pd
import requests

from .benchmark import prepare_cohort
from .data import load_measurements
from .surface import residue_table, summarize_surface

ROOT = Path(__file__).resolve().parents[2]
SITE = 'https://organizymedb.org'


def with_retry(operation):
    for attempt in range(3):
        try:
            return operation()
        except requests.RequestException as exc:
            if attempt==2:
                raise
            print('Network retry',attempt+1,type(exc).__name__,flush=True)
            time.sleep(2**attempt)


def structure_url(html):
    class Parser(HTMLParser):
        filename = None
        def handle_starttag(self, tag, attrs):
            if tag == 'body':
                self.filename = dict(attrs).get('data-pdb-file')
    parser = Parser()
    parser.feed(html)
    filename = parser.filename
    if not filename or '/' in filename or '\\' in filename or not filename.endswith('.cif'):
        raise ValueError('Missing or invalid published structure filename')
    return SITE+'/static/structures/'+quote(filename,safe='')


def candidate_ids(frame, catalog):
    lookup = {}
    for entry in catalog:
        key = (entry['enzyme_name'],entry['enzyme_species'])
        lookup.setdefault(key,set()).add(int(entry['protein_id']))
    return {sid:sorted({pid for row in group.itertuples()
                       for pid in lookup.get((row.enzyme_name,row.enzyme_species),set())})
            for sid,group in frame.groupby('sequence_id')}


def fetch_structure(protein_id, cache):
    # Some protein pages are tens of MB; stop after the body attribute, near the top.
    with requests.get(f'{SITE}/protein/{protein_id}/',stream=True,timeout=(10,25)) as response:
        response.raise_for_status()
        prefix = b''
        url = None
        for chunk in response.iter_content(4096):
            prefix += chunk
            try:
                url = structure_url(prefix.decode('utf-8',errors='replace'))
                break
            except ValueError:
                if len(prefix)>=131072:
                    raise ValueError('No structure attribute in first 128 KiB of page')
        if url is None:
            raise ValueError('No published structure link')
    path = cache/f'protein{protein_id}.cif'
    temporary = path.with_suffix('.cif.partial')
    with requests.get(url,stream=True,timeout=(10,25)) as response:
        response.raise_for_status()
        with temporary.open('wb') as handle:
            for chunk in response.iter_content(65536):
                handle.write(chunk)
    # Parse before promoting a complete response into the cache.
    from Bio.PDB import MMCIFParser
    MMCIFParser(QUIET=True).get_structure('download',str(temporary))
    temporary.replace(path)
    (cache/f'protein{protein_id}.source.json').write_text(json.dumps({
        'url':url,'protein_page':f'{SITE}/protein/{protein_id}/',
        'retrieved_utc':datetime.now(timezone.utc).isoformat()},indent=2))
    return path


def process_sequence(sequence_id, sequence, candidates, cache):
    record = dict(sequence_id=sequence_id,sequence_length=len(sequence),
                  sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                  candidate_protein_ids=candidates,status='unavailable',errors=[])
    for protein_id in candidates:
        try:
            path = cache/f'protein{protein_id}.cif'
            if not path.exists():
                path = fetch_structure(protein_id,cache)
            table,metadata = residue_table(path,sequence)
            assert ''.join(table.aa)==sequence
            assert len(table)==len(sequence)
            table.insert(0,'sequence_id',sequence_id)
            for cutoff in [10,20,30]:
                table[f'exposed_rsa{cutoff}'] = table.rsa.ge(cutoff/100)
            table['confident_plddt70'] = table.plddt.ge(70)
            features = dict(sequence_id=sequence_id,**summarize_surface(table))
            confident = table.loc[table.confident_plddt70]
            if not confident.empty and confident.rsa.ge(.2).any():
                features.update({'confident_'+key:value
                                 for key,value in summarize_surface(confident).items()})
            record.update(status='ok',protein_id=protein_id,structure_path=str(path),
                structure_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                mean_plddt=float(table.plddt.mean()),
                fraction_plddt_below70=float((~table.confident_plddt70).mean()),
                exposed_residues_rsa20=int(table.exposed_rsa20.sum()),**metadata)
            source = cache/f'protein{protein_id}.source.json'
            if source.exists():
                record['source'] = json.loads(source.read_text())
            return record,table,features
        except (ValueError,KeyError,AssertionError,OSError,requests.RequestException) as exc:
            record['errors'].append({'protein_id':protein_id,
                                     'error':f'{type(exc).__name__}: {exc}'})
    return record,None,None


def run(workers=2):
    out = ROOT/'reports/surface'
    cache = ROOT/'data/structures'
    out.mkdir(parents=True,exist_ok=True)
    cache.mkdir(parents=True,exist_ok=True)
    (out/'residues').mkdir(exist_ok=True)
    raw_path = ROOT/'data/raw/measurements.csv'
    frame,audit = prepare_cohort(load_measurements(raw_path))
    reference = json.loads((ROOT/'reports/family_benchmark/results.json').read_text())
    raw_hash = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    assert raw_hash==reference['audit']['raw_sha256']
    sequences = frame[['sequence_id','sequence']].drop_duplicates().sort_values('sequence_id')
    saved = pd.read_csv(ROOT/'reports/esm650m/source_sequences.csv')
    assert sequences.reset_index(drop=True).equals(saved)
    catalog_path = cache/'protein_catalog.json'
    if catalog_path.exists():
        catalog = json.loads(catalog_path.read_text())
    else:
        catalog = []
        page = 1
        while True:
            page_path = cache/f'protein_catalog_page{page}.json'
            if page_path.exists():
                data = json.loads(page_path.read_text())
            else:
                print('Fetching catalogue page',page,flush=True)
                def fetch_page():
                    response = requests.get(SITE+'/api/proteins/',
                        params={'page':page,'per_page':200},timeout=(10,25))
                    response.raise_for_status()
                    return response.json()
                data = with_retry(fetch_page)
                page_path.write_text(json.dumps(data,indent=2))
            catalog.extend(data['data'])
            if page>=data['pagination']['total_pages']:
                assert len(catalog)==data['pagination']['total']
                break
            page += 1
        assert len({r['protein_id'] for r in catalog})==len(catalog)
        catalog_path.write_text(json.dumps(catalog,indent=2))
    candidates = candidate_ids(frame,catalog)
    # Download each candidate once; keep network concurrency low.
    ids = sorted({pid for group in candidates.values() for pid in group})
    missing = [pid for pid in ids if not (cache/f'protein{pid}.cif').exists()]
    print('Cohort',len(sequences),'sequences; candidate structures',len(ids),
          '; uncached',len(missing),flush=True)
    download_errors = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(with_retry,lambda pid=pid:fetch_structure(pid,cache)):pid
                   for pid in missing}
        for i,future in enumerate(as_completed(futures),1):
            pid = futures[future]
            try:
                path = future.result()
            except Exception as exc:
                download_errors[str(pid)] = f'{type(exc).__name__}: {exc}'
            if i%10==0 or i==len(missing):
                print('Structures fetched',i,'/',len(missing),
                      'failures',len(download_errors),flush=True)
    records,features = [],[]
    for i,row in enumerate(sequences.itertuples(),1):
        available = [pid for pid in candidates[row.sequence_id]
                     if (cache/f'protein{pid}.cif').exists()]
        record,table,feature = process_sequence(row.sequence_id,row.sequence,available,cache)
        record['candidate_protein_ids'] = candidates[row.sequence_id]
        for pid in candidates[row.sequence_id]:
            if str(pid) in download_errors:
                record['errors'].append({'protein_id':pid,'error':download_errors[str(pid)]})
        records.append(record)
        if table is not None:
            table.to_csv(out/'residues'/f'{row.sequence_id}.csv',index=False)
            features.append(feature)
        (out/'manifest.json').write_text(json.dumps(records,indent=2))
        if i%25==0 or i==len(sequences):
            print('Surface calculated',i,'/',len(sequences),'matched',len(features),flush=True)
    pd.DataFrame(features).to_csv(out/'features.csv',index=False)
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        audit=audit,raw_sha256=raw_hash,requested=len(sequences),matched=len(features),
        unavailable=len(sequences)-len(features),download_errors=download_errors,
        exposed_definition='RSA >= 0.20; additional per-residue flags at 0.10 and 0.30',
        rsa_reference='FreeSASA default classifier reference areas; values above 1 retained',
        confidence='CA B-factor treated as pLDDT for database-provided predicted structures; >=70 sensitivity subset',
        geometry='All standard protein heavy atoms in the first provided model, without ligands or waters; not a verified biological assembly',
        versions={name:importlib.metadata.version(name) for name in ['biopython','freesasa','requests']})
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    print('SURFACE SUMMARY',json.dumps(summary),flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers',type=int,default=2)
    run(workers=parser.parse_args().workers)
