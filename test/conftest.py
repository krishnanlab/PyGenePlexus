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
    geneplexus.download.download_select_data(
        pytest.DATADIR,
        "All",
        "BioGRID",
        "Embedding",
        ["GO", "DisGeNet"],
        log_level="DEBUG",
    )
