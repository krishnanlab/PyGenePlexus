import os.path as osp
from typing import List
from typing import Tuple
from typing import Union
from urllib.parse import urljoin

import requests

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
from ._config.config import NET_SELECTION_TYPE
from ._config.config import NET_TYPE
from ._config.config import TASK_SELECTION_TYPE
from ._config.config import URL_AZURE


def download_select_data(
    data_dir: str,
    tasks: TASK_SELECTION_TYPE,
    networks: NET_SELECTION_TYPE = "All",
    features: FEATURE_SELECTION_TYPE = "All",
    GSCs: GSC_SELECTION_TYPE = "All",
):
    """Select subset of data to download.

    Args:
        data_dir (str): Location of data files.
        tasks: Task of interest ("IDConversion", "MachineLearning",
            "Similarities", "NetworkGraph"), accept multiple selection as a
            list. Do all the tasks if set to "All" (default: "All").
        networks: Networks of interest ("BioGRID", "STRING", "STRING-EXP",
            "GIANT-TN"), accept multiple selection as a list. Do all the
            networks if set to "All" (default: "All").
        features: Network features of interest ("Adjacency", "Embedding",
            "Influence"), accept multiple selection as a list. Do all the
            features if set to "All" (default: "All").
        GSCs: Gene set collection of interest ("GO", "DisGeNet"), accept
            multiple selection as a list. Do all the GSC if set to "All",
            (default: "All").

    """
    # Similarities and NetworkGraph will assume downloaded MachineLearning
    tasks, networks, features, GSCs = make_download_options_lists(tasks, networks, features, GSCs)
    all_files_to_do = []
    for atask in tasks:
        if atask == "IDconversion":
            all_files_to_do.extend(get_IDconversion_filenames())
        if atask == "MachineLearning":
            all_files_to_do.extend(get_MachineLearning_filenames(networks, features, GSCs))
        if atask == "Similarities":
            all_files_to_do.extend(get_Similarities_filenames(networks, features, GSCs))
        if atask == "NetworkGraph":
            all_files_to_do.extend(get_NetworkGraph_filenames(networks))

    all_files_to_do = list(set(all_files_to_do))
    download_from_azure(data_dir, all_files_to_do)


def download_from_azure(data_dir: str, files_to_do: List[str]):
    for afile in files_to_do:
        path = osp.join(data_dir, afile)
        if osp.exists(path):
            logger.info(f"File exists, skipping download: {path}")
        else:
            fn = urljoin(URL_AZURE, afile)
            logger.info(f"Downloading: {fn}")
            r = requests.get(fn)
            if r.ok:
                open(path, "wb").write(r.content)
            else:
                raise requests.exceptions.RequestException(r, fn)


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
    GSCs: GSC_SELECTION_TYPE,
) -> Tuple[List[str], List[NET_TYPE], List[FEATURE_TYPE], List[GSC_TYPE]]:
    args = (
        ("tasks", tasks, ALL_TASKS),
        ("network", networks, ALL_NETWORKS),
        ("feature", features, ALL_FEATURES),
        ("GSC", GSCs, ALL_GSCS),
    )
    return tuple(map(_make_download_options_list, *zip(*args)))  # type: ignore


def get_IDconversion_filenames() -> List[str]:
    files_to_do = []
    for line in util.get_all_filenames():
        if ("IDconversion" in line) or ("NodeOrder" in line):
            files_to_do.append(line)
    return files_to_do


def get_MachineLearning_filenames(
    networks: List[NET_TYPE],
    features: List[FEATURE_TYPE],
    GSCs: List[GSC_TYPE],
) -> List[str]:
    # TODO: switch GSCs and features position to make it consistent
    files_to_do = []
    for line in util.get_all_filenames():
        if "NodeOrder" in line:
            net_tmp = osp.splitext(line.split("Order_")[-1])[0]
            if net_tmp in networks:
                files_to_do.append(line)
        if ("universe.txt" in line) or ("GoodSets.json" in line):
            net_tmp = line.split("_")[2]
            GSC_tmp = line.split("_")[1]
            if (net_tmp in networks) and (GSC_tmp in GSCs):
                files_to_do.append(line)
        if "Data_" in line:
            feature_tmp = line.split("_")[1]
            net_tmp = osp.splitext(line.split("_")[2])[0]
            if (feature_tmp in features) and (net_tmp in networks):
                files_to_do.append(line)
        if ("Entrez-to-Name" in line) or ("Entrez-to-Symbol" in line):
            files_to_do.append(line)
    return files_to_do


def get_Similarities_filenames(
    networks: List[NET_TYPE],
    features: List[FEATURE_TYPE],
    GSCs: List[GSC_TYPE],
) -> List[str]:
    files_to_do = []
    for line in util.get_all_filenames():
        if "CorrectionMatrix_" in line:
            feature_tmp = osp.splitext(line.split("_")[-1])[0]
            net_tmp = line.split("_")[3]
            GSC_tmp = line.split("_")[1]
            if (net_tmp in networks) and (feature_tmp in features) and (GSC_tmp in GSCs):
                files_to_do.append(line)
        if "CorrectionMatrixOrder" in line:
            GSC_tmp = line.split("_")[1]
            net_tmp = osp.splitext(line.split("_")[2])[0]
            if (net_tmp in networks) and (GSC_tmp in GSCs):
                files_to_do.append(line)
        if "PreTrainedWeights" in line:
            net_tmp = line.split("_")[2]
            feature_tmp = osp.splitext(line.split("_")[-1])[0]
            if (net_tmp in networks) and (feature_tmp in features):
                files_to_do.append(line)
    return files_to_do


def get_NetworkGraph_filenames(networks: List[NET_TYPE]) -> List[str]:
    files_to_do = ["IDconversion_Homo-sapiens_Entrez-to-Symbol.json"]
    for line in util.get_all_filenames():
        if "Edgelist" in line:
            net_tmp = osp.splitext(line.split("_")[-1])[0]
            if net_tmp in networks:
                files_to_do.append(line)
    return files_to_do
