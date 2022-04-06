import os.path as osp
import pathlib

import pytest

import geneplexus


@pytest.fixture(scope="session")
def data(request):
    pytest.DATADIR = request.config.cache.makedir("download")
    pytest.HOMEDIR = pathlib.Path(__file__).parent.parent
    pytest.ANSWERDIR = osp.join(pytest.HOMEDIR, "test", "expected_result")
    pytest.GENELIST_PATH = osp.join(pytest.HOMEDIR, "example", "input_genes.txt")
    pytest.FILENAMES = [
        "CorrectionMatrix_DisGeNet_DisGeNet_BioGRID_Embedding.npy",
        "CorrectionMatrix_DisGeNet_GO_BioGRID_Embedding.npy",
        "CorrectionMatrix_GO_DisGeNet_BioGRID_Embedding.npy",
        "CorrectionMatrix_GO_GO_BioGRID_Embedding.npy",
        "CorrectionMatrixOrder_DisGeNet_BioGRID.txt",
        "CorrectionMatrixOrder_GO_BioGRID.txt",
        "Data_Embedding_BioGRID.npy",
        "Edgelist_BioGRID.edg",
        "GSC_DisGeNet_BioGRID_GoodSets.json",
        "GSC_DisGeNet_BioGRID_universe.txt",
        "GSC_GO_BioGRID_GoodSets.json",
        "GSC_GO_BioGRID_universe.txt",
        "IDconversion_Homo-sapiens_ENSG-to-Entrez.json",
        "IDconversion_Homo-sapiens_ENSP-to-Entrez.json",
        "IDconversion_Homo-sapiens_ENST-to-Entrez.json",
        "IDconversion_Homo-sapiens_Entrez-to-ENSG.json",
        "IDconversion_Homo-sapiens_Entrez-to-Name.json",
        "IDconversion_Homo-sapiens_Entrez-to-Symbol.json",
        "IDconversion_Homo-sapiens_Symbol-to-Entrez.json",
        "NodeOrder_BioGRID.txt",
        "NodeOrder_GIANT-TN.txt",
        "NodeOrder_STRING-EXP.txt",
        "NodeOrder_STRING.txt",
        "PreTrainedWeights_DisGeNet_BioGRID_Embedding.json",
        "PreTrainedWeights_GO_BioGRID_Embedding.json",
    ]

    geneplexus.download.download_select_data(
        pytest.DATADIR,
        "All",
        "BioGRID",
        "Embedding",
        ["GO", "DisGeNet"],
        log_level="DEBUG",
    )
