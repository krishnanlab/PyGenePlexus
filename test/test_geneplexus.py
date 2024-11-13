import numpy as np
import pytest

import geneplexus


# @pytest.mark.usefixtures("data")
@pytest.fixture(scope="module")
def gp():
    gp = geneplexus.GenePlexus(file_loc = pytest.DATADIR,
                               net_type = "STRING",
                               features = "SixSpeciesN2V",
                               sp_trn = "Human",
                               sp_res = "Mouse",
                               gsc_trn = "Combined",
                               gsc_res = "Combined")
    gp.load_genes(geneplexus.util.read_gene_list(pytest.GENELIST_PATH))
    return gp


@pytest.mark.parametrize("null_val", [None, -10])
@pytest.mark.parametrize("num_folds", [2, 3, 5])
@pytest.mark.parametrize("min_num_pos,cross_validate", [(100, True), (100, False), (200, True)])
# @pytest.mark.usefixtures("data")
def test_run_sl(gp, caplog, mocker, min_num_pos, num_folds, null_val, cross_validate):
    # Use random 5 dimensional vectors as features to speed up test
    mocker.patch(
        "geneplexus.util.load_gene_features",
        lambda w, x, y, z: np.random.random((30000, 5)),
    )

    gp.fit_and_predict(
        min_num_pos=min_num_pos,
        num_folds=num_folds,
        null_val=null_val,
        cross_validate=cross_validate,
    )

    if not cross_validate:
        assert "Skipping cross validation." in caplog.text
        assert gp.avgps == [null_val] * num_folds
    elif min_num_pos > len(gp.pos_genes_in_net):
        assert "Insufficient number of positive genes" in caplog.text
        assert f"{len(gp.pos_genes_in_net)} ({min_num_pos} needed)" in caplog.text
        assert gp.avgps == [null_val] * num_folds
    else:
        assert "Performing cross validation." in caplog.text

    assert len(gp.avgps) == num_folds
