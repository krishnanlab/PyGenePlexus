"""Data download module."""
import io
import os
import os.path as osp
import time
from concurrent.futures import ThreadPoolExecutor
from itertools import repeat
from threading import local
from typing import List
from typing import Tuple
from typing import Union
from urllib.parse import urljoin
from zipfile import ZipFile
import tarfile

import requests
from requests.sessions import Session

from . import util
from ._config import logger
from ._config.config import ALL_SPECIES
from ._config.config import SPECIES_SELECTION_TYPE
from ._config.config import SPECIES_TYPE
from ._config.config import LOG_LEVEL_TYPE
from ._config.config import MAX_RETRY
from ._config.config import URL_DICT
from ._config.logger_util import file_handler_context
from ._config.logger_util import stream_level_context
from .exception import DownloadError

thread_local = local()


def download_select_data(
    data_dir: str,
    species: SPECIES_SELECTION_TYPE = "All",
    data_loc: str = "Zenodo",
    n_jobs: int = 10,
    retry: bool = True,
    log_level: LOG_LEVEL_TYPE = "INFO",
):
    """Select species of data to download.

    Args:
        data_dir: Location of data files.
        species: Species of interest, accept multiple selection as a
            list. Do all the species if set to "All".
        n_jobs: Number of concurrent downloading threads.
        retry: If set to True, then retry downloading any missing file.

    """
    species = _get_species_list(species)
    with stream_level_context(logger, log_level):
        if data_loc == "Zenodo":
            logger.warn(
                f"Downloading data from Zenodo. This should take ~2 minutes per species"
                "but can vary greatly depending on download speeds",
            )
        for aspecies in species:
            if not _check_all_files(data_dir, aspecies):
                log_path = osp.join(data_dir, "download.log")
                logger.info(f"Start downloading data for {aspecies} and saving to: {data_dir}")
                fn_download = f"{aspecies}_data.tar.gz"
                with file_handler_context(logger, log_path, "DEBUG"):
                    _download_and_extract(data_dir, aspecies, fn_download, data_loc, retry)
                logger.info("Download completed.")


def _get_species_list(
    species: SPECIES_SELECTION_TYPE,
):
    if isinstance(species, str):
        if species == "All":
            species = ALL_SPECIES
        else:
            species = [species]
    elif  not isinstance(species, list):
        raise TypeError(f"Expected str type or list of str type, got {type(species)}")
    for i in species:
        if i not in ALL_SPECIES:
            raise ValueError(f"Unexpected species {i!r}")
    return species
    
def _check_all_files(
    data_dir: str,
    aspecies: str,
):
    fn_end = f"data_filenames_{aspecies}.txt"
    fn_full = osp.join(data_dir,fn_end)
    # check if filenames file is present
    if not osp.exists(fn_full):
        return False
    else:
        # if filenames file exsists, see if all files are present in data_dir
        with open(fn_full) as file:
            filenames = [line.rstrip() for line in file]
        files_found = [osp.basename(x) for x in os.listdir(data_dir)]
        files_missing = [x for x in filenames if x not in files_found]
        if len(files_missing) > 0:
            return False
        else:
            return True 

def _download_and_extract(data_dir, aspecies, fn_download, data_loc, retry):
    session = requests.Session()
    url = urljoin(URL_DICT[data_loc], fn_download)
    num_tries = 0
    while num_tries <= MAX_RETRY:
        num_tries += 1
        with session.get(url) as r:
            if r.ok:
                logger.debug(f"Response ok ({r!r}): {url=}")
                with tarfile.open(fileobj=io.BytesIO(r.content), mode="r:gz") as tf:
                    for member in tf.getmembers():
                        member.name = os.path.basename(member.name)
                        tf.extract(member,data_dir)    
                        logger.info(f"Downloaded {member.name}")
                if _check_all_files(data_dir, aspecies):    
                    break
                else:
                    logger.warning(f"Not all files downloaded, trying again")
            elif r.status_code == 429:  # Retry later
                t = r.headers["Retry-after"]
                logger.warning(f"Too many requests, waiting for {t} sec")
                time.sleep(int(t))
                continue
            else:
                raise requests.exceptions.RequestException(r, url)
        logger.critical("Session context closed, this should never happen!")
    else:
        raise DownloadError(f"Failed to download from {url} ({MAX_RETRY=})")
