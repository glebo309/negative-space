"""Conservative data preparation and statistical checks for transfer benchmarks."""
import re
import numpy as np
import pandas as pd

from .data import select_primary_endpoint


def connected_groups(nodes, edges):
    parent = {node: node for node in nodes}

    def root(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for left, right in edges:
        left, right = root(left), root(right)
        if left != right:
            parent[max(left, right)] = min(left, right)
    return {node: root(node) for node in nodes}


def assignment_table(sequences, groupings):
    exported = sequences.copy()
    for split, mapping in groupings.items():
        column = 'sequence_group' if split == 'sequence' else split
        exported[column] = exported.sequence_id.map(mapping)
        if exported[column].isna().any():
            raise ValueError('Incomplete grouping assignments')
    return exported


def parse_condition(value, stage):
    text = '' if pd.isna(value) else str(value)
    match = re.search(rf'{stage}:\s*([^,;]+)', text, re.I)
    if not match:
        return None
    number = re.fullmatch(r'\s*([-+]?\d+(?:\.\d+)?)\s*(?:°\s*C)?\s*', match[1])
    return float(number[1]) if number else None


def incubation_hours(value):
    text = '' if pd.isna(value) else str(value)
    match = re.search(r'incubation\s*\(\s*(\d+(?:\.\d+)?)\s*(seconds?|minutes?|hours?|days?)\b', text, re.I)
    if not match:
        return None
    factors = {'second': 1/3600, 'minute': 1/60, 'hour': 1, 'day': 24}
    return float(match[1]) * factors[match[2].lower().rstrip('s')]


def prepare_cohort(raw):
    frame = select_primary_endpoint(raw)
    audit = {'initial_primary_rows': len(frame)}
    frame = frame.loc[frame.sequence.str.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+')].copy()
    audit['standard_sequence_rows'] = len(frame)
    control = pd.to_numeric(frame.aqueous_control, errors='coerce')
    frame = frame.loc[control.eq(100) & frame.measured_value.ge(0) &
                      frame.solvent_volume.gt(0) & frame.solvent_volume.le(100) &
                      frame.doi.notna() & frame.solvent_name.notna()].copy()
    audit['control100_valid_rows'] = len(frame)
    # Remove only identical records, except the database row identifier.
    subset = [c for c in frame.columns if c != 'measurement_id']
    dedup = frame.drop_duplicates(subset=subset).copy()
    audit['duplicate_rows_removed'] = len(frame) - len(dedup)
    frame = dedup.reset_index(drop=True)
    for column, suffix in [('temperature','temperature_c'), ('ph','ph')]:
        for stage in ['Incubation','Assay']:
            frame[f'{stage.lower()}_{suffix}'] = frame[column].map(lambda x: parse_condition(x,stage))
    frame['incubation_hours'] = frame.procedure.map(incubation_hours)
    frame['solvent_logp'] = pd.to_numeric(frame.solvent_logp,errors='coerce')
    frame['target_log1p'] = np.log1p(frame.measured_value)
    sequences = sorted(frame.sequence.unique())
    ids = {seq: f'seq{i:04d}' for i, seq in enumerate(sequences)}
    frame['sequence_id'] = frame.sequence.map(ids)
    audit.update(rows=len(frame),sequences=len(sequences),solvents=frame.solvent_name.nunique(),publications=frame.doi.nunique())
    return frame, audit


def paired_cluster_ci(scores, n_boot=5000, seed=11):
    if scores.empty:
        return {'n_sequences':0, 'n_clusters':0, 'delta':None, 'ci95':None}
    values = scores.assign(delta=scores.composition-scores.baseline)
    groups = values.groupby('cluster').delta.agg(['sum','count'])
    rng = np.random.default_rng(seed)
    draws = rng.integers(0,len(groups),size=(n_boot,len(groups)))
    differences = groups['sum'].to_numpy()[draws].sum(axis=1)/groups['count'].to_numpy()[draws].sum(axis=1)
    return {'n_sequences':len(scores), 'n_clusters':len(groups),
            'delta':float(values.delta.mean()), 'ci95':np.quantile(differences,[.025,.975]).tolist()}


def rank_scores(frame, matched=False):
    rows = []
    keys = ['sequence_id']
    if matched:
        keys += [c for c in ['doi','solvent_volume','procedure','temperature','ph',
            'assay_solution','substrates','substrate_concentrations','cofactor','shaking','comments'] if c in frame]
    for _, group in frame.groupby(keys,dropna=False):
        ranked = group
        if matched:
            ranked = group.groupby('solvent_name')[['measured_value','baseline','composition']].mean()
        if len(ranked)<3 or ranked.measured_value.nunique()<2:
            continue
        # A constant prediction has no ranking information; score it zero.
        correlations = {name: (float(ranked.measured_value.corr(ranked[name],method='spearman'))
            if ranked[name].nunique()>1 else 0.) for name in ['baseline','composition']}
        rows.append(dict(sequence_id=group.sequence_id.iloc[0],cluster=group.cluster.iloc[0], **correlations))
    if not rows:
        return pd.DataFrame(columns=['sequence_id','cluster','baseline','composition','n_series'])
    scores = pd.DataFrame(rows)
    return scores.groupby(['sequence_id','cluster'],as_index=False).agg(
        baseline=('baseline','mean'),composition=('composition','mean'),n_series=('baseline','size'))
