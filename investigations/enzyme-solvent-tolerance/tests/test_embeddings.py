import numpy as np
import pandas as pd
import pytest
import torch

from enzsolv.embeddings import residue_mean, sequence_chunks, join_embeddings


def test_residue_mean_excludes_padding_and_special_tokens():
    hidden = torch.tensor([[[99.,99.],[1.,3.],[3.,5.],[88.,88.],[77.,77.]]])
    attention = torch.tensor([[1,1,1,1,0]])
    special = torch.tensor([[1,0,0,1,1]])
    assert torch.allclose(residue_mean(hidden,attention,special),torch.tensor([[2.,4.]]))


def test_chunks_cover_long_sequence_without_truncation():
    sequence = 'A'*1022+'C'*80
    chunks = sequence_chunks(sequence)
    assert list(map(len,chunks)) == [1022,80]
    assert ''.join(chunks)==sequence
    with pytest.raises(ValueError):
        sequence_chunks('')


def test_join_embeddings_uses_ids_and_rejects_missing_records():
    frame = pd.DataFrame({'sequence_id':['b','a','b']})
    result = join_embeddings(frame,['a','b'],np.array([[1.,2.],[3.,4.]]))
    assert result.esm_000.tolist()==[3.,1.,3.]
    with pytest.raises(ValueError,match='Missing'):
        join_embeddings(frame,['a'],np.array([[1.,2.]]))
    with pytest.raises(ValueError,match='Duplicate'):
        join_embeddings(frame,['a','a'],np.array([[1.,2.],[3.,4.]]))
