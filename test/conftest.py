import os.path as osp
import pathlib

import pytest

import geneplexus


@pytest.fixture(scope="session")
def data(request):
    pytest.DATADIR = request.config.cache.makedir("download")
    pytest.RESULTSDIR = request.config.cache.makedir("results")
    pytest.CLIRESULTSDIR = request.config.cache.makedir("cli_results")
    pytest.HOMEDIR = pathlib.Path(__file__).parent.parent
    pytest.ANSWERDIR = osp.join(pytest.HOMEDIR, "test", "expected_result")
    pytest.GENELIST_PATH = osp.join(pytest.HOMEDIR, "example", "input_genes.txt")
    pytest.FILENAMES = [
        "Data__Human__SixSpeciesN2V__STRING.npy",
        "IDconversion__Human__ENSG-to-Entrez.json",
        "IDconversion__Human__ENSP-to-Entrez.json",
        "IDconversion__Human__ENST-to-Entrez.json",
        "IDconversion__Human__Entrez-to-Name.json",
        "IDconversion__Human__Entrez-to-Symbol.json",
        "IDconversion__Human__Symbol-to-Entrez.json",
        "NodeOrder__Human__BioGRID.txt",
        "NodeOrder__Human__IMP.txt",
        "NodeOrder__Human__STRING.txt",
        "BioMart__Human__Mouse.json",
        "Data__Mouse__SixSpeciesN2V__STRING.npy",
        "Edgelist__Mouse__STRING.edg",
        "GSC__Human__Combined__STRING.json",
        "GSC__Mouse__Combined__STRING.json",
        "IDconversion__Mouse__Entrez-to-Name.json",
        "IDconversion__Mouse__Entrez-to-Symbol.json",
        "NodeOrder__Mouse__STRING.txt",
        "PreTrainedWeights__Mouse__Combined__STRING__SixSpeciesN2V.json",
    ]

    geneplexus.download.download_pytest_data(
        pytest.DATADIR,
        log_level="DEBUG",
    )
