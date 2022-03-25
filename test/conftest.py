import pytest

import geneplexus


@pytest.fixture(scope="session")
def data(request):
    pytest.DATADIR = request.config.cache.makedir("download")
    geneplexus.download.download_select_data(
        pytest.DATADIR,
        "All",
        "BioGRID",
        "Embedding",
        ["GO", "DisGeNet"],
    )
