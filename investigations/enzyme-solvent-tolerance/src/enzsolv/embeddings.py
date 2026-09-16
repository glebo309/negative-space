"""Frozen ESM representations with explicit residue pooling and ID alignment."""
import numpy as np
import pandas as pd
import torch


def residue_mean(hidden, attention, special):
    mask = (attention.bool() & ~special.bool()).unsqueeze(-1)
    count = mask.sum(dim=1)
    if (count == 0).any():
        raise ValueError('No residues to pool')
    return (hidden * mask).sum(dim=1) / count


def sequence_chunks(sequence, limit=1022):
    if not sequence or limit < 1:
        raise ValueError('A nonempty sequence and positive limit are required')
    return [sequence[i:i+limit] for i in range(0,len(sequence),limit)]


def join_embeddings(frame, ids, vectors):
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate embedding identifiers')
    lookup = pd.DataFrame(vectors,index=ids,columns=[f'esm_{i:03d}' for i in range(vectors.shape[1])])
    if not set(frame.sequence_id).issubset(lookup.index):
        raise ValueError('Missing embedding identifiers')
    return frame.join(lookup,on='sequence_id',validate='many_to_one')
