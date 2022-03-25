import pytest

import geneplexus


@pytest.fixture(scope="class")
def data(request):
    pytest.datadir = request.config.cache.makedir("download")
    geneplexus.download.download_select_data(
        pytest.datadir,
        "All",
        "BioGRID",
        "Embedding",
        ["GO", "DisGeNet"],
    )
