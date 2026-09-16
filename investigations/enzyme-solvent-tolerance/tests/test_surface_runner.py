import numpy as np
import pandas as pd
import pytest

from enzsolv.run_surface_benchmark import fit_oof


def tiny_rows():
    return pd.DataFrame({'measurement_id':range(6), 'sequence_id':['a']*3+['b']*3,
        'cluster':['a']*3+['b']*3, 'doi':['paper1']*3+['paper2']*3,
        'fold':[0]*3+[1]*3, 'solvent_name':['A','B','C']*2,
        'solvent_volume':[5.,10.,20.]*2, 'target_log1p':np.log1p([10,20,30,15,25,35])})


def test_oof_predicts_every_row_and_cache_is_keyed_by_inputs(tmp_path):
    rows=tiny_rows()
    predicted=fit_oof(rows,['solvent_volume'],tmp_path/'test.npz',trees=5)
    assert predicted.shape==(6,)
    assert np.isfinite(predicted).all()
    assert np.array_equal(predicted,fit_oof(rows,['solvent_volume'],tmp_path/'test.npz',trees=5))
    changed=rows.copy()
    changed.loc[0,'solvent_volume']=99.
    with pytest.raises(ValueError,match='Cache'):
        fit_oof(changed,['solvent_volume'],tmp_path/'test.npz',trees=5)


def test_oof_rejects_cluster_leakage_and_publication_leakage(tmp_path):
    rows=tiny_rows()
    rows.loc[3:,'cluster']='a'
    with pytest.raises(ValueError,match='cluster'):
        fit_oof(rows,['solvent_volume'],tmp_path/'test.npz',trees=5)
    rows=tiny_rows()
    rows.loc[3:,'doi']='paper1'
    with pytest.raises(ValueError,match='publication'):
        fit_oof(rows,['solvent_volume'],tmp_path/'other.npz',trees=5,publication_holdout=True)
