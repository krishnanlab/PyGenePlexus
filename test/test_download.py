import shutil
import tempfile
from urllib.parse import urljoin

import pytest

from geneplexus import download
from geneplexus._config.config import URL_DATA
from geneplexus.exception import DownloadError


@pytest.mark.usefixtures("data")
def test_retries(requests_mock, tmpdir):
    for filename in pytest.FILENAMES:
        requests_mock.get(urljoin(URL_DATA, f"{filename}.zip"), status_code=429)

    with pytest.raises(DownloadError) as excinfo:
        download._download_from_url(tmpdir, pytest.FILENAMES)
    assert str(excinfo.value) == "Failed to download all required files (MAX_RETRY=10)"
