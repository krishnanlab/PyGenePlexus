import os.path as osp
import pathlib

import pytest

import geneplexus


@pytest.fixture(scope="session")
def data(request):
    pytest.DATADIR = request.config.cache.makedir("download")
    pytest.HOMEDIR = pathlib.Path(__file__).parent.parent
    pytest.ANSWERDIR = osp.join(pytest.HOMEDIR, "test", "expected_result")
    geneplexus.download.download_select_data(
        pytest.DATADIR,
        "All",
        "BioGRID",
        "Embedding",
        ["GO", "DisGeNet"],
    )
