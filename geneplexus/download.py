import os.path as osp
from urllib.parse import urljoin

import requests

from . import util
from ._config import config
from ._config import logger


def download_select_data(fp_data, tasks="All", networks="All", features="All", GSCs="All"):
    # Similarities and NetworkGraph will assume downloaded MachineLearning
    tasks, networks, features, GSCs = make_download_options_lists(tasks, networks, features, GSCs)
    all_files_to_do = []
    for atask in tasks:
        if atask == "IDconversion":
            files_to_do = get_IDconversion_filenames()
            all_files_to_do = all_files_to_do + files_to_do
        if atask == "MachineLearning":
            files_to_do = get_MachineLearning_filenames(networks, GSCs, features)
            all_files_to_do = all_files_to_do + files_to_do
        if atask == "Similarities":
            files_to_do = get_Similarities_filenames(networks, features, GSCs)
            all_files_to_do = all_files_to_do + files_to_do
        if atask == "NetworkGraph":
            files_to_do = get_NetworkGraph_filenames(networks)
            all_files_to_do = all_files_to_do + files_to_do

    all_files_to_do = list(set(all_files_to_do))
    download_from_azure(fp_data, all_files_to_do)


def download_from_azure(fp_data, files_to_do):
    for afile in files_to_do:
        path = osp.join(fp_data, afile)
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


def make_download_options_lists(tasks, networks, features, GSCs):
    if isinstance(tasks, str):
        if tasks == "All":
            tasks = config.ALL_TASKS
        elif tasks in config.ALL_TASKS:
            tasks = [tasks]
        else:
            raise ValueError(f"Unexpected task: {tasks!r}")
    if isinstance(networks, str):
        if networks == "All":
            networks = config.ALL_NETWORKS
        elif networks in config.ALL_NETWORKS:
            networks = [networks]
        else:
            raise ValueError(f"Unexpected network: {tasks!r}")
    if isinstance(features, str):
        if features == "All":
            features = config.ALL_FEATURES
        elif features in config.ALL_FEATURES:
            features = [features]
        else:
            raise ValueError(f"Unexpected feature: {features!r}")
    if isinstance(GSCs, str):
        if GSCs == "All":
            GSCs = config.ALL_GSCS
        elif GSCs in config.ALL_GSCS:
            GSCs = [GSCs]
        else:
            raise ValueError(f"Unexpected GSC: {GSCs!r}")

    return tasks, networks, features, GSCs


def get_IDconversion_filenames():
    files_to_do = []
    for line in util.get_filenames():
        if ("IDconversion" in line) or ("NodeOrder" in line):
            files_to_do.append(line)
    return files_to_do


def get_MachineLearning_filenames(networks, GSCs, features):
    files_to_do = []
    for line in util.get_filenames():
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
    for line in util.get_filenames():
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
    for line in util.get_filenames():
        if "Edgelist" in line:
            net_tmp = osp.splitext(line.split("_")[-1])[0]
            if net_tmp in networks:
                files_to_do.append(line)
    return files_to_do
