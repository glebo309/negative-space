import numpy as np
import pandas as pd
import pytest

from enzsolv.surface import summarize_surface, residue_table


def test_surface_features_use_exposed_residues_and_area():
    residues = pd.DataFrame({'aa':list('ADKF'), 'rsa':[.8,.4,.05,.01],
        'asa':[80.,40.,5.,1.], 'plddt':[90.]*4})
    features = summarize_surface(residues)
    assert features['surface_exposed_fraction']==.5
    assert features['surface_aa_A']==.5
    assert features['surface_aa_D']==.5
    assert features['surface_aa_K']==0
    assert np.isclose(features['surface_area_acidic_fraction'],40/126)


def test_surface_chain_requires_exact_sequence_identity(tmp_path):
    # A synthetic residue with a complete alanine heavy-atom set.
    from Bio.PDB import Atom, Chain, Model, Residue, Structure
    from Bio.PDB.mmcifio import MMCIFIO
    structure = Structure.Structure('x')
    model = Model.Model(0)
    chain = Chain.Chain('A')
    residue = Residue.Residue((' ',1,' '),'ALA','')
    for i,(name,coord,element) in enumerate([
        ('N',[0,0,0],'N'),('CA',[1.45,0,0],'C'),('C',[2,1.4,0],'C'),
        ('O',[2,2.4,0],'O'),('CB',[1.5,-.5,1.4],'C')]):
        residue.add(Atom.Atom(name,np.array(coord,dtype=float),90,1,' ',name,i,element=element))
    chain.add(residue); model.add(chain); structure.add(model)
    path=tmp_path/'a.cif'
    writer=MMCIFIO(); writer.set_structure(structure); writer.save(str(path))
    table,metadata = residue_table(path,'A')
    assert table.aa.tolist()==['A']
    assert table.asa.iloc[0]>0
    assert table.rsa.iloc[0]>0
    assert metadata['matched_chain']=='A'
    with pytest.raises(ValueError,match='exact'):
        residue_table(path,'C')
