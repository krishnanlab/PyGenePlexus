"""Data download module."""
import io
import os
import os.path as osp
import shutil
import tarfile
import time
from typing import List
from typing import Tuple
from typing import Union
from urllib.parse import urljoin

import pystow
import requests

from . import util
from ._config import logger
from ._config.config import ALL_SPECIES
from ._config.config import LOG_LEVEL_TYPE
from ._config.config import MAX_RETRY
from ._config.config import SPECIES_SELECTION_TYPE
from ._config.config import SPECIES_TYPE
from ._config.config import URL_DICT
from ._config.logger_util import file_handler_context
from ._config.logger_util import stream_level_context
from .exception import DownloadError


def download_select_data(
    file_loc: str = None,
    species: SPECIES_SELECTION_TYPE = "All",
    data_loc: str = "ZenodoAPI",
    num_retries: int = MAX_RETRY,
    log_level: LOG_LEVEL_TYPE = "INFO",
):
    """Select species of data to download.

    Args:
        file_loc: Location to save data files to . if not specified, set to default
            data path ``~/.data/geneplexus``
        species: Species of interest, accept multiple selection as a
            list. Do all the species if set to "All".
        data_loc: the remote system where to look for the data
        num_retries: Number of times to retry downloading a file.
        log_level: Level to set the logger

    """
    if file_loc is None:
        file_loc = str(pystow.join("geneplexus"))
    else:
        file_loc = util.normexpand(file_loc)
    species = _get_species_list(species)
    with stream_level_context(logger, log_level):
        for aspecies in species:
            if not _check_all_files(file_loc, aspecies):
                if data_loc in ["Zenodo", "ZenodoAPI"]:
                    logger.warning(
                        f"Downloading {aspecies} data from Zenodo. This should take ~2 "
                        "minutes per species but can vary greatly depending on download speeds. "
                        "If Zenodo download is hanging for > 5 minutes per attempt, it might be best "
                        "to stop and restart the PyGenePlexus download function.",
                    )
                log_path = osp.join(file_loc, "download.log")
                logger.info(f"Start downloading data for {aspecies} and saving to: {file_loc}")
                fn_download = f"{aspecies}_data.tar.gz"
                if data_loc == "ZenodoAPI":
                    fn_download = f"{fn_download}/content"
                with file_handler_context(logger, log_path, "DEBUG"):
                    _download_and_extract(file_loc, aspecies, fn_download, data_loc, num_retries)
                logger.info("Download completed.")
            else:
                logger.warning(
                    f"Files already downloaded for {aspecies}",
                )


def _get_species_list(
    species: SPECIES_SELECTION_TYPE,
):
    if isinstance(species, str):
        if species == "All":
            species = ALL_SPECIES
        else:
            species = [species]
    elif not isinstance(species, list):
        raise TypeError(f"Expected str type or list of str type, got {type(species)}")
    for i in species:
        if i not in ALL_SPECIES:
            raise ValueError(f"Unexpected species {i!r}")
    return species


def _check_all_files(
    file_loc: str,
    file_cat: str,
):
    fn_end = f"data_filenames_{file_cat}.txt"
    fn_full = osp.join(file_loc, fn_end)
    # check if filenames file is present
    if not osp.exists(fn_full):
        return False
    else:
        # if filenames file exsists, see if all files are present in file_loc
        with open(fn_full) as file:
            filenames = [line.rstrip() for line in file]
        files_found = [osp.basename(x) for x in os.listdir(file_loc)]
        files_missing = [x for x in filenames if x not in files_found]
        if len(files_missing) > 0:
            return False
        else:
            return True


def _download_and_extract(file_loc, file_cat, fn_download, data_loc, num_retries):
    url = urljoin(URL_DICT[data_loc], fn_download)
    num_tries = 0
    while num_tries <= num_retries - 1:
        num_tries += 1
        logger.info(f"On attempt {num_tries} of {num_retries} for downloading the data")
        try:
            with requests.get(url, timeout=2) as r:
                if r.ok:
                    logger.debug(f"Response ok ({r!r}): {url=}")
                    with tarfile.open(fileobj=io.BytesIO(r.content), mode="r:gz") as tf:
                        for member in tf.getmembers():
                            member.name = os.path.basename(member.name)
                            tf.extract(member, file_loc)
                            logger.info(f"Downloaded {member.name}")
                    try:
                        shutil.rmtree(osp.join(file_loc, f"{file_cat}_data"))
                    except FileNotFoundError:
                        pass
                    if _check_all_files(file_loc, file_cat):
                        break
                    else:
                        logger.warning(f"Not all files downloaded, trying again")
                elif r.status_code == 429:  # Retry later
                    t = r.headers["Retry-after"]
                    logger.warning(f"Too many requests, waiting for {t} sec")
                    time.sleep(int(t))
                    continue
                else:
                    logger.info("An unknown error occured")
                    continue
        except:
            logger.info("An error occured during download (probably a connection timeout)")
            continue
        logger.critical("Session context closed, this should never happen!")
    else:
        raise DownloadError(f"Failed to download from {url} ({num_retries=})")


def download_pytest_data(
    file_loc: str = None,
    data_loc: str = "ZenodoAPI",
    num_retries: int = MAX_RETRY,
    log_level: LOG_LEVEL_TYPE = "INFO",
):
    """Download data for pytests.

    Args:
        file_loc: Location of data files.
        data_loc: the remote system where to look for the data
        num_retries: Number of times to retry downloading a file.
        log_level: Level to set the logger

    """
    if file_loc is None:
        file_loc = str(pystow.join("geneplexus"))
    else:
        file_loc = util.normexpand(file_loc)
    with stream_level_context(logger, log_level):
        if not _check_all_files(file_loc, "pytest"):
            if data_loc in ["Zenodo", "ZenodoAPI"]:
                logger.warning(
                    f"Downloading pytest data from Zenodo. This should take ~2 "
                    "minutes but can vary greatly depending on download speeds. "
                    "If Zenodo download is hanging for > 5 minutes per attempt, it might be best "
                    "to stop and restart the PyGenePlexus download function.",
                )
            log_path = osp.join(file_loc, "download.log")
            logger.info(f"Start downloading pytest data and saving to: {file_loc}")
            fn_download = "pytest_data.tar.gz"
            if data_loc == "ZenodoAPI":
                fn_download = f"{fn_download}/content"
            with file_handler_context(logger, log_path, "DEBUG"):
                _download_and_extract(file_loc, "pytest", fn_download, data_loc, num_retries)
            logger.info("Download completed.")
        else:
            logger.warning(
                f"Files already downloaded for pytest",
            )
