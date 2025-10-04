import os
import os.path as osp
import pathlib
import time

import numpy as np

import geneplexus

# The first step is for a user to supply a gene list
input_genes = [
    "ARL6",
    "BBS1",
    "BBS10",
    "BBS12",
    "BBS2",
    "BBS4",
    "BBS5",
    "BBS7",
    "BBS9",
    "CCDC28B",
    "CEP290",
    "KIF7",
    "MKKS",
    "MKS1",
    "TRIM32",
    "TTC8",
    "WDPCP",
]

# bardet biedl syndrome gene set
input_genes = ["6457","7037","57403","3134","50807","93343","11311","8766","5584","137492","998","30011","5337","3312","155",
               "10015","55738","57132","153","116986","163","11267","1950","3559","6714","84249","2066","29924","1213","30846",
               "84612","440073","2060","3303","3561","9101","51160","56904","3304","23527","5878","3560","7189","3949","92421",
               "26286","5979","9922","11031","116983","2261","9230","5867","64145","867","57154","84313","3577","116987","10617",
               "1436","200576","83737","23396","3310","5590","3133","382","6456","30845","868","2264","5868","84440","116984",
               "5869","23624","22841","161","23096","5338","652614","84552","51028","55616","9829","3815","29082","9135","23362",
               "9146","128866","156","8218","89853","154","64744","9525","84364","9727","23550","8853","1956","8395","6455",
               "64411","5156","51100","8027","408","3305","51534","2868","9744","3106","51652","3265","27243","10938","60682",
               "157","26056","10059","2321","80230","1173","1175","160","3306","3135","1234","2149","8411","3791","51510","23327",
               "409","11059","3579","27183","8396","1601","1211","3480","9815","26119","64750","26052","4914","25978","8394","1212",
               "30844","131890","79720","7251","50855","116985","5662","2870","10193","1785","155382","652799","22905","3105",
               "55048","10254","55040","7852","1759","4193","2869","2065","6011","4734","28964","4233","80223","79643","3107",
               "2263","56288"]


# Set up directories
homedir = pathlib.Path(__file__).absolute().parent
datadir = osp.join(homedir, "data")
outdir = osp.join(homedir, "result")
os.makedirs(datadir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)

# Get the data from URL
geneplexus.download.download_select_data(
    datadir,
    species=["Human", "Mouse"],
)

# Run through the pipeline
# First initialize the geneplexus object
gp = geneplexus.GenePlexus(
    file_loc=datadir,
    net_type="STRING",
    features="SixSpeciesN2V",
    sp_trn="Human",
    sp_res=["Human", "Mouse"],
    gsc_trn="Combined",
    gsc_res=["Combined", "Combined"],
)

# Load the input genes into the class and set up positives/negatives
gp.load_genes(input_genes)

# Run the steps of the pipeline, including the option clustering step
gp.cluster_input()
gp.fit()
gp.predict()
gp.make_sim_dfs()
gp.make_small_edgelist()

# get human gene prediction results for full input gene set model
print(gp.model_info["All-Genes"].results["Human-Combined"].df_probs)
# get similarties of trainied model to other models trained with human annotations
print(gp.model_info["All-Genes"].results["Human-Combined"].df_sim)
# get network connections for the top 50 human genes predcited using full input gene set model
print(gp.model_info["All-Genes"].results["Human-Combined"].df_edge_sym)

# get mouse gene prediction results for cluster 1 gene set model
print(gp.model_info["Cluster-01"].results["Mouse-Combined"].df_probs)
# get similarties of trainied model to other models trained with mouse annotations
print(gp.model_info["Cluster-01"].results["Mouse-Combined"].df_sim)
# get network connections for the top 50 mouse genes predcited using cluster 1 gene set model
print(gp.model_info["Cluster-01"].results["Mouse-Combined"].df_edge_sym)

# get log2(auPRC/prior) metric for the full input gene set model
print(gp.model_info["All-Genes"].avgps)

# save the class. If output_dir=None will try to save to ~/.data/geneplexus_outputs/results
gp.save_class(output_dir = outdir)