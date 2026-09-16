import numpy as np
import pandas as pd

from enzsolv.benchmark import connected_groups, parse_condition, incubation_hours, prepare_cohort, paired_cluster_ci, rank_scores, assignment_table


def test_connected_groups_preserves_transitive_similarity():
    groups = connected_groups(['a', 'b', 'c', 'd'], [('a', 'b'), ('b', 'c')])
    assert groups['a'] == groups['c']
    assert groups['d'] != groups['a']


def test_condition_parser_never_turns_ranges_or_units_into_point_estimates():
    assert parse_condition('Incubation: 25°C, Assay: 37°C', 'Incubation') == 25
    assert parse_condition('Incubation: ?, Assay: 8', 'Incubation') is None
    assert parse_condition('Incubation: 20-30°C, Assay: 37°C', 'Incubation') is None
    assert parse_condition('30°C', 'Incubation') is None
    assert parse_condition('Incubation: Room temperature, Assay: 30°C', 'Incubation') is None


def test_incubation_time_requires_unambiguous_duration():
    assert incubation_hours('Activity after incubation (30 minutes at 50°C), assay 10 minutes') == .5
    assert incubation_hours('Activity after incubation (2 hours at 25°C)') == 2
    assert incubation_hours('Activity after incubation (1 day at 25°C)') == 24
    assert incubation_hours('Activity after incubation (1-2 hours at 25°C)') is None
    assert incubation_hours('Activity after incubation (overnight at 25°C)') is None


def test_cohort_preserves_source_sequences_and_different_procedures():
    row = dict(doi='d',enzyme_name='e',enzyme_species='s',mutation='WT',sequence='AACC',
        property='Stability - Incubation',solvent_name='Methanol',solvent_volume='10',
        measured_value='80',units='%',aqueous_control='100',ph='Incubation: 7, Assay: 8',
        temperature='Incubation: 25°C, Assay: 37°C',procedure='Activity after incubation (1 hour at 25°C)',
        solvent_logp='1')
    data = pd.DataFrame([row, row, {**row,'procedure':'Activity after incubation (2 hours at 25°C)'},
        {**row,'aqueous_control':'0'}, {**row,'sequence':'AA CX'}, {**row,'aqueous_control':'50'}])
    selected, audit = prepare_cohort(data)
    assert len(selected) == 2
    assert audit['duplicate_rows_removed'] == 1
    assert data.loc[4,'sequence'] == 'AA CX'
    assert set(selected.sequence) == {'AACC'}
    assert set(selected.incubation_hours) == {1,2}


def test_bootstrap_pairs_models_and_resamples_clusters():
    scores = pd.DataFrame({'cluster':['a','a','b'], 'baseline':[.1,.2,.3], 'composition':[.3,.4,.5]})
    result = paired_cluster_ci(scores, n_boot=100)
    assert np.isclose(result['delta'], .2)
    assert np.allclose(result['ci95'], [.2,.2])


def test_solvent_ranking_does_not_mix_concentrations():
    frame = pd.DataFrame(dict(sequence_id=['a']*6, cluster=['g']*6, doi=['d']*6,
        solvent_volume=[10]*3+[20]*3, solvent_name=['x','y','z']*2,
        measured_value=[1,2,3,4,5,6], baseline=[3,2,1,6,5,4], composition=[1,2,3,4,5,6]))
    scores = rank_scores(frame, matched=True)
    assert len(scores)==1
    assert np.isclose(scores.iloc[0].baseline,-1)
    assert np.isclose(scores.iloc[0].composition,1)
    assert scores.iloc[0].n_series==2
    frame['solvent_volume'] = range(6)
    assert rank_scores(frame, matched=True).empty


def test_assignment_table_cannot_overwrite_amino_acid_sequences():
    sequences = pd.DataFrame({'sequence_id':['a','b'], 'sequence':['AACC','GGTT']})
    exported = assignment_table(sequences,{'sequence':{'a':'a','b':'b'},'identity30':{'a':'a','b':'a'}})
    assert exported.sequence.tolist() == ['AACC','GGTT']
    assert exported.sequence_group.tolist() == ['a','b']
    assert exported.identity30.tolist() == ['a','a']
