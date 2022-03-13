import os.path as osp
from typing import List
from typing import Tuple
from typing import Union
from urllib.parse import urljoin

import requests

from . import config
from . import util
from ._config import logger


def download_select_data(
    data_dir: str,
    tasks: Union[str, List[str]] = "All",
    networks: Union[str, List[str]] = "All",
    features: Union[str, List[str]] = "All",
    GSCs: Union[str, List[str]] = "All",
):
    # Similarities and NetworkGraph will assume downloaded MachineLearning
    tasks, networks, features, GSCs = make_download_options_lists(tasks, networks, features, GSCs)
    all_files_to_do = []
    for atask in tasks:
        if atask == "IDconversion":
            all_files_to_do.extend(get_IDconversion_filenames())
        if atask == "MachineLearning":
            all_files_to_do.extend(get_MachineLearning_filenames(networks, GSCs, features))
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
            fn = urljoin(config.URL_AZURE, afile)
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
    tasks: Union[str, List[str]],
    networks: Union[str, List[str]],
    features: Union[str, List[str]],
    GSCs: Union[str, List[str]],
) -> Tuple[List[str], List[str], List[str], List[str]]:
    return map(
        _make_download_options_list,
        *zip(
            ("tasks", tasks, config.ALL_TASKS),
            ("network", networks, config.ALL_NETWORKS),
            ("feature", features, config.ALL_FEATURES),
            ("GSC", GSCs, config.ALL_GSCS),
        ),
    )


def get_IDconversion_filenames():
    files_to_do = []
    for line in util.get_all_filenames():
        if ("IDconversion" in line) or ("NodeOrder" in line):
            files_to_do.append(line)
    return files_to_do


def get_MachineLearning_filenames(networks, GSCs, features):
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


def get_Similarities_filenames(networks, features, GSCs):
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


def get_NetworkGraph_filenames(networks):
    files_to_do = ["IDconversion_Homo-sapiens_Entrez-to-Symbol.json"]
    for line in util.get_all_filenames():
        if "Edgelist" in line:
            net_tmp = osp.splitext(line.split("_")[-1])[0]
            if net_tmp in networks:
                files_to_do.append(line)
    return files_to_do
