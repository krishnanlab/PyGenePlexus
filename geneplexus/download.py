"""Data download module."""
import io
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

import requests
from requests.sessions import Session

from . import util
from ._config import logger
from ._config.config import ALL_FEATURES
from ._config.config import ALL_GSCS
from ._config.config import ALL_NETWORKS
from ._config.config import ALL_TASKS
from ._config.config import FEATURE_SELECTION_TYPE
from ._config.config import FEATURE_TYPE
from ._config.config import GSC_SELECTION_TYPE
from ._config.config import GSC_TYPE
from ._config.config import LOG_LEVEL_TYPE
from ._config.config import MAX_RETRY
from ._config.config import NET_SELECTION_TYPE
from ._config.config import NET_TYPE
from ._config.config import TASK_SELECTION_TYPE
from ._config.config import TASK_TYPE
from ._config.config import URL_DATA
from ._config.logger_util import log_level_context
from .exception import DownloadError

thread_local = local()


def download_select_data(
    data_dir: str,
    tasks: TASK_SELECTION_TYPE = "All",
    networks: NET_SELECTION_TYPE = "All",
    features: FEATURE_SELECTION_TYPE = "All",
    gscs: GSC_SELECTION_TYPE = "All",
    n_jobs: int = 10,
    retry: bool = True,
    log_level: LOG_LEVEL_TYPE = "INFO",
):
    """Select subset of data to download.

    Args:
        data_dir: Location of data files.
        tasks: Task of interest, accept multiple selection as a list. Do all
            the tasks if set to "All".
        networks: Networks of interest, accept multiple selection as a list. Do
            all the networks if set to "All".
        features: Network features of interest, accept multiple selection as a
            list. Do all the features if set to "All".
        gscs: Gene set collection of interest, accept multiple selection as a
            list. Do all the GSC if set to "All".
        n_jobs: Number of concurrent downloading threads.
        retry: If set to True, then retry downloading any missing file.

    """
    # Similarities and NetworkGraph will assume downloaded MachineLearning
    tasks, networks, features, gscs = make_download_options_lists(tasks, networks, features, gscs)
    all_files_to_do = []
    for atask in tasks:
        if atask == "IDconversion":
            all_files_to_do.extend(get_id_conversion_filenames())
        if atask == "MachineLearning":
            all_files_to_do.extend(get_machine_learning_filenames(networks, features, gscs))
        if atask == "Similarities":
            all_files_to_do.extend(get_similarities_filenames(networks, features, gscs))
        if atask == "NetworkGraph":
            all_files_to_do.extend(get_network_filenames(networks))
        if atask == "OriginalGSCs":
            all_files_to_do.extend(get_original_gscs_filenames())

    with log_level_context(logger, log_level):
        if features != ["Embedding"]:
            logger.warn(
                f"Downloading data type {features!r} may take a while (~10min "
                "to an hour depending on the downloadspeed)",
            )
        files_to_download = _get_files_to_download(data_dir, list(set(all_files_to_do)))
        if len(files_to_download) > 0:
            logger.info(f"Total number of files to download: {len(files_to_download)}")
            logger.info(f"Start downloading data and saving to: {data_dir}")
            _download_from_url(data_dir, files_to_download, n_jobs, retry)
            logger.info(f"Download completed.")


def _get_session() -> Session:
    if not hasattr(thread_local, "session"):
        thread_local.session = requests.Session()
        logger.debug(f"Acquired thread local session {thread_local.session!r}")
    return thread_local.session


def _download_file(file: str, data_dir: str):
    session = _get_session()
    url = urljoin(URL_DATA, f"{file}.zip")
    logger.debug(f"Thread started: {url=}, {session=}")
    num_tries = 1
    while num_tries <= MAX_RETRY:
        with session.get(url) as r:
            if r.ok:
                logger.debug(f"Response ok ({r!r}): {url=}")
                ZipFile(io.BytesIO(r.content)).extractall(data_dir)
                break
            elif r.status_code == 429:  # Retry later
                t = r.headers["Retry-after"]
                logger.warning(f"Too many requests, waiting for {t} sec")
                time.sleep(int(t))
                num_tries += 1
                continue
            else:
                raise requests.exceptions.RequestException(r, url)
        logger.critical("Session context closed, this should never happen!")
    else:
        raise DownloadError(f"Failed to download from {url} ({MAX_RETRY=})")
    logger.info(f"Downloaded {file}")


def _get_files_to_download(
    data_dir: str,
    files: List[str],
    silent: bool = False,
) -> List[str]:
    files_to_download = []
    for file in files:
        path = osp.join(data_dir, file)
        if osp.exists(path):
            if not silent:
                logger.debug(f"File exists, skipping download: {path}")
        else:
            files_to_download.append(file)
    return files_to_download


def _download_from_url(
    data_dir: str,
    files_to_do: List[str],
    n_jobs: int = 10,
    retry: bool = True,
    retry_count: int = 0,
):
    """Download file using the base url.

    Args:
        data_dir: Location of data files.
        files_to_do: List of files to download from the the url.
        n_jobs: Number of concurrent downloading threads.
        retry: If set to True, then retry downloading any missing file.
        retry_count: (DO NOT MODIFY) Counting the number of retries.

    """
    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        executor.map(_download_file, files_to_do, repeat(data_dir))

    missed_files = _get_files_to_download(data_dir, files_to_do, silent=True)
    if missed_files and retry:
        if retry_count >= MAX_RETRY:
            raise DownloadError(f"Failed to download all required files ({MAX_RETRY=})")
        missed = "".join(f"\n\t{i}" for i in missed_files)
        logger.warning(
            f"Failed to download the following files, retrying...{missed}",
        )
        _download_from_url(data_dir, missed_files, n_jobs, True, retry_count + 1)


def _make_download_options_list(
    name: str,
    opts: Union[str, List[str]],
    check_list: List[str],
) -> List[str]:
    if isinstance(opts, str):
        opts = check_list if opts == "All" else [opts]
    elif not isinstance(opts, list):
        raise TypeError(f"Expcted str type or list of str type, got {type(opts)}")

    for i in opts:
        if i not in check_list:
            raise ValueError(f"Unexpected {name}: {opts!r}")

    return opts


def make_download_options_lists(
    tasks: TASK_SELECTION_TYPE,
    networks: NET_SELECTION_TYPE,
    features: FEATURE_SELECTION_TYPE,
    gscs: GSC_SELECTION_TYPE,
) -> Tuple[List[TASK_TYPE], List[NET_TYPE], List[FEATURE_TYPE], List[GSC_TYPE]]:
    """Compile a list of files to download based on the selections."""
    args = (
        ("tasks", tasks, ALL_TASKS),
        ("network", networks, ALL_NETWORKS),
        ("feature", features, ALL_FEATURES),
        ("gsc", gscs, ALL_GSCS),
    )
    return tuple(map(_make_download_options_list, *zip(*args)))  # type: ignore


def get_id_conversion_filenames() -> List[str]:
    """Get gene ID conversion file names."""
    files_to_do = []
    for line in util.get_all_filenames():
        if ("IDconversion" in line) or ("NodeOrder" in line):
            files_to_do.append(line)
    return files_to_do


def get_machine_learning_filenames(
    networks: List[NET_TYPE],
    features: List[FEATURE_TYPE],
    gscs: List[GSC_TYPE],
) -> List[str]:
    """Get dataset file names."""
    files_to_do = []
    for line in util.get_all_filenames():
        if "NodeOrder" in line:
            net_tmp = osp.splitext(line.split("Order_")[-1])[0]
            if net_tmp in networks:
                files_to_do.append(line)
        if ("universe.txt" in line) or ("GoodSets.json" in line):
            net_tmp = line.split("_")[2]
            gsc_tmp = line.split("_")[1]
            if (net_tmp in networks) and (gsc_tmp in gscs):
                files_to_do.append(line)
        if "Data_" in line:
            feature_tmp = line.split("_")[1]
            net_tmp = osp.splitext(line.split("_")[2])[0]
            if (feature_tmp in features) and (net_tmp in networks):
                files_to_do.append(line)
        if ("Entrez-to-Name" in line) or ("Entrez-to-Symbol" in line):
            files_to_do.append(line)
    return files_to_do


def get_similarities_filenames(
    networks: List[NET_TYPE],
    features: List[FEATURE_TYPE],
    gscs: List[GSC_TYPE],
) -> List[str]:
    """Get pretrained model similarity file names."""
    files_to_do = []
    for line in util.get_all_filenames():
        if "CorrectionMatrix_" in line:
            feature_tmp = osp.splitext(line.split("_")[-1])[0]
            net_tmp = line.split("_")[3]
            gsc_tmp = line.split("_")[1]
            if (net_tmp in networks) and (feature_tmp in features) and (gsc_tmp in gscs):
                files_to_do.append(line)
        if "CorrectionMatrixOrder" in line:
            gsc_tmp = line.split("_")[1]
            net_tmp = osp.splitext(line.split("_")[2])[0]
            if (net_tmp in networks) and (gsc_tmp in gscs):
                files_to_do.append(line)
        if "PreTrainedWeights" in line:
            net_tmp = line.split("_")[2]
            feature_tmp = osp.splitext(line.split("_")[-1])[0]
            if (net_tmp in networks) and (feature_tmp in features):
                files_to_do.append(line)
    return files_to_do


def get_network_filenames(networks: List[NET_TYPE]) -> List[str]:
    """Get network file names."""
    files_to_do = ["IDconversion_Homo-sapiens_Entrez-to-Symbol.json"]
    for line in util.get_all_filenames():
        if "Edgelist" in line:
            net_tmp = osp.splitext(line.split("_")[-1])[0]
            if net_tmp in networks:
                files_to_do.append(line)
    return files_to_do


def get_original_gscs_filenames() -> List[str]:
    """Get original GSC file names."""
    files_to_do = []
    for line in util.get_all_filenames():
        if "GSCOriginal" in line:
            files_to_do.append(line)
    return files_to_do
