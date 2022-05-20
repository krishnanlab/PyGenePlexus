import pytest

import geneplexus


@pytest.mark.usefixtures("data")
@pytest.fixture(scope="module")
def gp():
    gp = geneplexus.GenePlexus(pytest.DATADIR, "BioGRID", "Embedding", "GO")
    gp.load_genes(geneplexus.util.read_gene_list(pytest.GENELIST_PATH))
    return gp


@pytest.mark.parametrize("null_val", [-10, -20])
@pytest.mark.parametrize("num_folds", [2, 3, 5])
@pytest.mark.parametrize("min_num_pos,cross_validate", [(100, True), (100, False), (200, True)])
@pytest.mark.usefixtures("data")
def test_run_sl(gp, caplog, min_num_pos, num_folds, null_val, cross_validate):
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
