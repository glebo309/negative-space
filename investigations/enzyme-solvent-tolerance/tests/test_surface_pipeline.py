import pandas as pd
import pytest

from enzsolv.run_surface import candidate_ids, structure_url, process_sequence


def test_network_retry_recovers_timeout_and_does_not_retry_bad_data(monkeypatch):
    import requests
    import time
    import enzsolv.run_surface as module
    monkeypatch.setattr(time,'sleep',lambda seconds:None)
    calls = []
    def operation():
        calls.append(1)
        if len(calls)<3:
            raise requests.Timeout('temporary interruption')
        return 'ok'
    assert module.with_retry(operation)=='ok'
    assert len(calls)==3
    with pytest.raises(ValueError,match='bad structure'):
        module.with_retry(lambda: (_ for _ in ()).throw(ValueError('bad structure')))


def test_structure_url_uses_published_attribute_and_escapes_filename():
    html = '<body data-pdb-file="Lipase A_(gene_lipA).cif">'
    assert structure_url(html) == 'https://organizymedb.org/static/structures/Lipase%20A_%28gene_lipA%29.cif'
    for html in ['<body>', '<body data-pdb-file="../secret.cif">',
                 '<body data-pdb-file="https://elsewhere/a.cif">']:
        with pytest.raises(ValueError):
            structure_url(html)


def test_candidate_ids_require_exact_name_and_species_and_cover_all_sequences():
    frame = pd.DataFrame({'sequence_id':['s1','s1','s2'], 'enzyme_name':['X','X','Y'],
                          'enzyme_species':['A','A','B']})
    catalog = [{'protein_id':7,'enzyme_name':'X','enzyme_species':'A'},
               {'protein_id':8,'enzyme_name':'X','enzyme_species':'Other'}]
    assert candidate_ids(frame,catalog) == {'s1':[7], 's2':[]}


def test_process_sequence_records_missing_structures_without_fake_features(tmp_path):
    manifest, table, features = process_sequence('s1','ACDE',[],tmp_path)
    assert manifest['status']=='unavailable'
    assert manifest['sequence_id']=='s1'
    assert table is None and features is None


def test_process_sequence_uses_cached_cif_and_rejects_mismatch(tmp_path,monkeypatch):
    import enzsolv.run_surface as module
    path = tmp_path/'protein7.cif'
    path.write_text('cached structure')
    def fetch(*args,**kwargs):
        raise AssertionError('Cached structures must not be downloaded again')
    monkeypatch.setattr(module,'fetch_structure',fetch)
    def calculate(path,sequence):
        if sequence!='AD':
            raise ValueError('No exact sequence match in structure')
        return pd.DataFrame({'position':[1,2], 'aa':['A','D'], 'asa':[80.,40.],
             'rsa':[.8,.4], 'plddt':[90.,50.]}), {'matched_chain':'A'}
    monkeypatch.setattr(module,'residue_table',calculate)
    manifest, table, features = process_sequence('s1','AD',[7],tmp_path)
    assert manifest['status']=='ok'
    assert table.exposed_rsa20.tolist()==[True,True]
    assert table.confident_plddt70.tolist()==[True,False]
    assert features['surface_aa_D']==.5
    assert features['confident_surface_aa_D']==0
    assert manifest['sequence_length']==2
    manifest, table, features = process_sequence('s2','AC',[7],tmp_path)
    assert manifest['status']=='unavailable'
    assert 'exact sequence' in manifest['errors'][0]['error']
    assert table is None and features is None
