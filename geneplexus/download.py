import os.path as osp

import requests


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
            print(f"The following file already exsists so skipping download: {path}")
        else:
            FN_Azure = f"https://mancusogeneplexusstorage.blob.core.windows.net/mancusoblob2/{afile}"
            print(f"Downloading the follwing file: {FN_Azure}")
            r = requests.get(FN_Azure)
            open(path, "wb").write(r.content)


def make_download_options_lists(tasks, networks, features, GSCs):
    all_tasks = ["IDconversion", "MachineLearning", "Similarities", "NetworkGraph"]
    all_networks = ["BioGRID", "STRING", "STRING-EXP", "GIANT-TN"]
    all_features = ["Adjacency", "Embedding", "Influence"]
    all_GSCs = ["GO", "DisGeNet"]

    if isinstance(tasks, str):
        if tasks == "All":
            tasks = all_tasks
        elif tasks in all_tasks:
            tasks = [tasks]
        else:
            raise ValueError(f"Unexpected task: {tasks!r}")
    if isinstance(networks, str):
        if networks == "All":
            networks = all_networks
        elif networks in all_networks:
            networks = [networks]
        else:
            raise ValueError(f"Unexpected network: {tasks!r}")
    if isinstance(features, str):
        if features == "All":
            features = all_features
        elif features in all_features:
            features = [features]
        else:
            raise ValueError(f"Unexpected feature: {features!r}")
    if isinstance(GSCs, str):
        if GSCs == "All":
            GSCs = all_GSCs
        elif GSCs in all_GSCs:
            GSCs = [GSCs]
        else:
            raise ValueError(f"Unexpected GSC: {GSCs!r}")

    return tasks, networks, features, GSCs


def get_IDconversion_filenames():
    files_to_do = []
    with open("data_filenames.txt", "r") as f:
        for line in f:
            line = line.strip()
            if ("IDconversion" in line) or ("NodeOrder" in line):
                files_to_do.append(line)
    return files_to_do


def get_MachineLearning_filenames(networks, GSCs, features):
    files_to_do = []
    with open("data_filenames.txt", "r") as f:
        for line in f:
            line = line.strip()
            if "NodeOrder" in line:
                net_tmp = line.split("Order_")[-1].split(".tx")[0]
                if net_tmp in networks:
                    files_to_do.append(line)
            if ("universe.txt" in line) or ("GoodSets.pickle" in line):
                net_tmp = line.split("_")[2]
                GSC_tmp = line.split("_")[1]
                if (net_tmp in networks) and (GSC_tmp in GSCs):
                    files_to_do.append(line)
            if "Data_" in line:
                feature_tmp = line.split("_")[1]
                net_tmp = line.split("_")[2].split(".n")[0]
                if (feature_tmp in features) and (net_tmp in networks):
                    files_to_do.append(line)
            if ("Entrez-to-Name" in line) or ("Entrez-to-Symbol" in line):
                files_to_do.append(line)
    return files_to_do


def get_Similarities_filenames(networks, features, GSCs):
    files_to_do = []
    with open("data_filenames.txt", "r") as f:
        for line in f:
            line = line.strip()
            if "CorrectionMatrix_" in line:
                feature_tmp = line.split("_")[-1].split(".n")[0]
                net_tmp = line.split("_")[3]
                GSC_tmp = line.split("_")[1]
                if (net_tmp in networks) and (feature_tmp in features) and (GSC_tmp in GSCs):
                    files_to_do.append(line)
            if "CorrectionMatrixOrder" in line:
                GSC_tmp = line.split("_")[1]
                net_tmp = line.split("_")[2].split(".t")[0]
                if (net_tmp in networks) and (GSC_tmp in GSCs):
                    files_to_do.append(line)
            if "PreTrainedWeights" in line:
                net_tmp = line.split("_")[2]
                feature_tmp = line.split("_")[-1].split(".pic")[0]
                if (net_tmp in networks) and (feature_tmp in features):
                    files_to_do.append(line)
    return files_to_do


def get_NetworkGraph_filenames(networks):
    files_to_do = ["IDconversion_Homo-sapiens_Entrez-to-Symbol.pickle"]
    with open("data_filenames.txt", "r") as f:
        for line in f:
            line = line.strip()
            if "Edgelist" in line:
                net_tmp = line.split("_")[-1].split(".ed")[0]
                if net_tmp in networks:
                    files_to_do.append(line)
    return files_to_do
