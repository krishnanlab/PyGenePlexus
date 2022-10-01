from urllib.parse import urljoin

import pytest

from geneplexus import download
from geneplexus._config import logger
from geneplexus._config.config import URL_DICT
from geneplexus._config.logger_util import log_level_context
from geneplexus.exception import DownloadError


@pytest.mark.usefixtures("data")
def test_retries(requests_mock, tmpdir):
    for filename in pytest.FILENAMES:
        requests_mock.get(urljoin(URL_DICT["Zenodo"], f"{filename}.zip"), status_code=429)

    with log_level_context(logger, "CRITICAL"):
        with pytest.raises(DownloadError) as excinfo:
            download._download_from_url(tmpdir, pytest.FILENAMES, "Zenodo")
    assert str(excinfo.value) == "Failed to download all required files (MAX_RETRY=10)"
