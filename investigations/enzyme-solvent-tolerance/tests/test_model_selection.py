from enzsolv.run_esm import model_spec


def test_model_selection_keeps_small_and_large_results_separate():
    small = model_spec('35m')
    large = model_spec('650m')
    assert small[0] != large[0]
    assert small[0].name=='esm35m'
    assert large[0].name=='esm650m'
    assert small[1]=='facebook/esm2_t12_35M_UR50D'
    assert large[1]=='facebook/esm2_t33_650M_UR50D'
    assert len(small[2])==40
    assert len(large[2])==40
