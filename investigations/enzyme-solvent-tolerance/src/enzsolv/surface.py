"""Calculate residue solvent accessibility from sequence-verified structures."""
import string
import numpy as np
import pandas as pd
import freesasa
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def residue_table(path, expected_sequence):
    structure = MMCIFParser(QUIET=True).get_structure('enzyme',str(path))[0]
    chains = []
    for chain in structure:
        residues = [r for r in chain if r.id[0]==' ' and seq1(r.resname) in AMINO_ACIDS]
        if residues:
            chains.append((chain.id,residues,''.join(seq1(r.resname) for r in residues)))
    matches = [i for i,(_,_,sequence) in enumerate(chains) if sequence==expected_sequence]
    if not matches:
        raise ValueError('No exact sequence match in structure')
    labels = string.ascii_uppercase+string.ascii_lowercase+string.digits
    if len(chains)>len(labels):
        raise ValueError('Too many chains for FreeSASA labels')
    geometry = freesasa.Structure()
    for i,(_,residues,_) in enumerate(chains):
        for position,residue in enumerate(residues,1):
            for atom in residue:
                if atom.element not in {'H','D'}:
                    geometry.addAtom(f' {atom.name:<3}'[:4],residue.resname,str(position),labels[i],
                        *[float(x) for x in atom.coord])
    result = freesasa.calc(geometry,freesasa.Parameters({'probe-radius':1.4,'n-threads':1}))
    matched = matches[0]
    areas = result.residueAreas()[labels[matched]]
    rows = []
    for position,residue in enumerate(chains[matched][1],1):
        area = areas[str(position)]
        if not area.hasRelativeAreas or 'CA' not in residue:
            raise ValueError('Missing reference area or alpha carbon')
        rows.append(dict(position=position,aa=seq1(residue.resname),asa=area.total,
            rsa=area.relativeTotal,polar_asa=area.polar,apolar_asa=area.apolar,
            plddt=float(residue['CA'].bfactor)))
    return pd.DataFrame(rows),dict(matched_chain=chains[matched][0],
        matching_chains=len(matches),protein_chains=len(chains),probe_radius_angstrom=1.4,
        geometry='First model, all standard protein heavy atoms; solvent and ligands excluded')


def summarize_surface(residues, threshold=.2):
    exposed = residues.rsa.ge(threshold)
    if not exposed.any() or not np.isfinite(residues[['asa','rsa']]).all().all():
        raise ValueError('Invalid or empty exposed surface')
    features = {'surface_exposed_fraction':float(exposed.mean()),
        'surface_mean_rsa':float(residues.rsa.mean()),
        'surface_asa_per_residue':float(residues.asa.sum()/len(residues))}
    for aa in AMINO_ACIDS:
        features[f'surface_aa_{aa}'] = float(residues.loc[exposed,'aa'].eq(aa).mean())
    for name,letters in {'acidic':'DE','basic':'KR','hydrophobic':'AVILMFWY','polar':'STNQ'}.items():
        mask = residues.aa.isin(list(letters))
        features[f'surface_area_{name}_fraction'] = float(residues.loc[mask,'asa'].sum()/residues.asa.sum())
    return features
