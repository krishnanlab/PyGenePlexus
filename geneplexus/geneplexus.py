import os
import pickle
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine
from scipy.stats import hypergeom
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

import utls


class GenePlexus:
    def __init__(
        self,
        file_loc="/Users/christophermancuso/Documents/DataSets/Geneplexus_data/",
        network="BioGRID",
        features="Embedding",
        GSC="GO",
    ):
        self.file_loc = file_loc
        self.network = network
        self.features = features
        self.GSC = GSC

    def load_genes(self, input_genes):
        input_genes = [item.upper() for item in input_genes]
        self.input_genes = input_genes

    def convert_to_Entrez(self):
        convert_IDs, df_convert_out = utls.intial_ID_convert(self.input_genes, self.file_loc)
        df_convert_out, table_summary, input_count = utls.make_validation_df(df_convert_out, self.file_loc)
        self.convert_IDs = convert_IDs
        self.df_convert_out = df_convert_out
        self.table_summary = table_summary
        self.input_count = input_count
        return self.df_convert_out

    def set_params(self, net_type, features, GSC):
        self.net_type = net_type
        self.features = features
        self.GSC = GSC

    def get_pos_and_neg_genes(self):
        pos_genes_in_net, genes_not_in_net, net_genes = utls.get_genes_in_network(
            self.file_loc,
            self.net_type,
            self.convert_IDs,
        )
        self.pos_genes_in_net = pos_genes_in_net
        self.genes_not_in_net = genes_not_in_net
        self.net_genes = net_genes
        negative_genes = utls.get_negatives(self.file_loc, self.net_type, self.GSC, self.pos_genes_in_net)
        self.negative_genes = negative_genes
        return self.pos_genes_in_net, self.negative_genes, self.net_genes

    def fit_and_predict(self):
        mdl_weights, probs, avgps = utls.run_SL(
            self.file_loc,
            self.net_type,
            self.features,
            self.pos_genes_in_net,
            self.negative_genes,
            self.net_genes,
        )
        self.mdl_weights = mdl_weights
        self.probs = probs
        self.avgps = avgps
        df_probs = utls.make_prob_df(
            self.file_loc,
            self.net_genes,
            self.probs,
            self.pos_genes_in_net,
            self.negative_genes,
        )
        self.df_probs = df_probs
        return self.mdl_weights, self.df_probs, self.avgps

    def make_sim_dfs(self):
        df_sim_GO, df_sim_Dis, weights_GO, weights_Dis = utls.make_sim_dfs(
            self.file_loc,
            self.mdl_weights,
            self.GSC,
            self.net_type,
            self.features,
        )
        self.df_sim_GO = df_sim_GO
        self.df_sim_Dis = df_sim_Dis
        self.weights_GO = weights_GO
        self.weights_Dis = weights_Dis
        return self.df_sim_GO, self.df_sim_Dis, self.weights_GO, self.weights_Dis

    def make_small_edgelist(self, num_nodes=50):
        df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = utls.make_small_edgelist(
            self.file_loc,
            self.df_probs,
            self.net_type,
            num_nodes=50,
        )
        self.df_edge = df_edge
        self.isolated_genes = isolated_genes
        self.df_edge_sym = df_edge_sym
        self.isolated_genes_sym = isolated_genes_sym
        return self.df_edge, self.isolated_genes, self.df_edge_sym, self.isolated_genes_sym

    def alter_validation_df(self):
        df_convert_out_subset, positive_genes = utls.alter_validation_df(
            self.df_convert_out,
            self.table_summary,
            self.net_type,
        )
        self.df_convert_out_subset = df_convert_out_subset
        self.positive_genes = positive_genes
        return self.df_convert_out_subset, self.positive_genes


def download_IDconversion_data(fp_data):
    files_to_do = utls.get_IDconversion_filenames()
    utls.download_from_azure(fp_data, files_to_do)


def download_all_data(fp_data):
    with open("data_filenames.txt", "r") as f:
        for line in f:
            line = line.strip()
            if os.path.exists(fp_data + line):
                print("The following file already exsists so skipping download", fp_data + line)
            else:
                FN_Azure = "https://mancusogeneplexusstorage.blob.core.windows.net/mancusoblob2/%s" % line
                print("Downloading file from", FN_Azure)
                r = requests.get(FN_Azure)
                open(fp_data + line, "wb").write(r.content)


def download_select_data(fp_data, tasks="All", networks="All", features="All", GSCs="All"):
    # Similarities and NetworkGraph will assume downloaded MachineLearning
    tasks, networks, features, GSCs = utls.make_download_options_lists(tasks, networks, features, GSCs)
    all_files_to_do = []
    for atask in tasks:
        if atask == "IDconversion":
            files_to_do = utls.get_IDconversion_filenames()
            all_files_to_do = all_files_to_do + files_to_do
        if atask == "MachineLearning":
            files_to_do = utls.get_MachineLearning_filenames(networks, GSCs, features)
            all_files_to_do = all_files_to_do + files_to_do
        if atask == "Similarities":
            files_to_do = utls.get_Similarities_filenames(networks, features, GSCs)
            all_files_to_do = all_files_to_do + files_to_do
        if atask == "NetworkGraph":
            files_to_do = utls.get_NetworkGraph_filenames(networks)
            all_files_to_do = all_files_to_do + files_to_do

    all_files_to_do = list(set(all_files_to_do))
    print("The number of files to download is", len(files_to_do))
    # utls.download_from_azure(fp_data,all_files_to_do)
