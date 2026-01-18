import os.path as osp
import shutil

import numpy as np
import pytest

import geneplexus
from geneplexus.exception import NoPositivesError


# @pytest.mark.usefixtures("data")
@pytest.fixture(scope="module")
def gp():
    gp = geneplexus.GenePlexus(
        file_loc=pytest.DATADIR,
        net_type="STRING",
        features="SixSpeciesN2V",
        sp_trn="Human",
        sp_res="Mouse",
        gsc_trn="Combined",
        gsc_res="Combined",
    )
    gp.load_genes(geneplexus.util.read_gene_list(pytest.GENELIST_PATH))
    return gp


@pytest.mark.parametrize("null_val", [None, -10])
@pytest.mark.parametrize("num_folds", [2, 3, 5])
@pytest.mark.parametrize("min_num_pos_cv,cross_validate", [(100, True), (100, False), (200, True)])
@pytest.mark.parametrize(
    "min_num_pos,excepted_error_message",
    [(10, None), (200, "There were not enough positive genes to train the model")],
)  # current example geneset has 183 genes
def test_run_sl(
    gp,
    caplog,
    mocker,
    min_num_pos,
    min_num_pos_cv,
    num_folds,
    null_val,
    cross_validate,
    excepted_error_message,
):
    # Use random 5 dimensional vectors as features to speed up test
    mocker.patch(
        "geneplexus.util.load_gene_features",
        lambda w, x, y, z: np.random.random((30000, 5)),
    )

    try:
        gp.fit(
            min_num_pos=min_num_pos,
            min_num_pos_cv=min_num_pos_cv,
            num_folds=num_folds,
            null_val=null_val,
            cross_validate=cross_validate,
        )

        if not cross_validate:
            assert "Skipping cross validation." in caplog.text
            assert gp.model_info["All-Genes"].avgps == [null_val] * num_folds
        elif min_num_pos_cv > len(gp.model_info["All-Genes"].pos_genes_in_net):
            assert "Insufficient number of positive genes" in caplog.text
            assert f"{len(gp.model_info['All-Genes'].pos_genes_in_net)} ({min_num_pos_cv} needed)" in caplog.text
            assert gp.model_info["All-Genes"].avgps == [null_val] * num_folds
        else:
            assert "Performing cross validation." in caplog.text

        assert len(gp.model_info["All-Genes"].avgps) == num_folds

    except NoPositivesError as e:
        assert excepted_error_message in str(e)


@pytest.mark.parametrize("sp_res", ["Human", ["Human", "Mouse"]])
def test_res_params(
    caplog,
    mocker,
    sp_res,
):
    # Use random 5 dimensional vectors as features to speed up test
    mocker.patch(
        "geneplexus.util.load_gene_features",
        lambda w, x, y, z: np.random.random((30000, 5)),
    )

    gp2 = geneplexus.GenePlexus(
        file_loc=pytest.DATADIR,
        net_type="STRING",
        features="SixSpeciesN2V",
        sp_trn="Human",
        sp_res=sp_res,
        gsc_trn="Combined",
        gsc_res="Combined",
    )
    gp2.load_genes(geneplexus.util.read_gene_list(pytest.GENELIST_PATH))
    gp2.fit()
    gp2.predict()

    if sp_res == "Human":
        assert "Human-Combined" in gp2.model_info["All-Genes"].results
    elif sp_res == ["Human", "Mouse"]:
        assert "Human-Combined" in gp2.model_info["All-Genes"].results
        assert "Mouse-Combined" in gp2.model_info["All-Genes"].results


@pytest.mark.parametrize("clust_method", ["louvain"])
def test_clustering(
    gp,
    caplog,
    mocker,
    clust_method,
):
    # Use random 5 dimensional vectors as features to speed up test
    mocker.patch(
        "geneplexus.util.load_gene_features",
        lambda w, x, y, z: np.random.random((30000, 5)),
    )
    extra_el = osp.join(pytest.HOMEDIR, "test", "extra_test_data", "Edgelist__Human__STRING.edg")
    shutil.copy(extra_el, pytest.DATADIR)

    gp.cluster_input(clust_method=clust_method)
    gp.fit()
    gp.predict()

    assert "Cluster-01" in gp.model_info


if __name__ == "__main__":
    unittest.main()
