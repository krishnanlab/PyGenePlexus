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

# Set up directories
homedir = pathlib.Path(__file__).absolute().parent
datadir = osp.join(homedir, "data")
outdir = osp.join(homedir, "result")
os.makedirs(datadir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)

# Get the data from URL
geneplexus.download.download_select_data(
    datadir,
    tasks="All",
    networks="STRING",
    features="Embedding",
    gscs=["GO", "DisGeNet"],
)

# Run through the pipeline
# First initialize the geneplexus object
myclass = geneplexus.GenePlexus(datadir, "STRING", "Embedding", "DisGeNet")

# Load the input genes into the class and set up positives/negatives
myclass.load_genes(input_genes)

# This trains the model and predcits on every gene in the network
mdl_weights, df_probs, avgps = myclass.fit_and_predict()

# The makes the tables that have the model weight similarity to other models
# trained on known GO and DisGeNet sets
df_sim_GO, df_sim_Dis, weights_GO, weights_Dis = myclass.make_sim_dfs()

# Return an edgelist
df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = myclass.make_small_edgelist(num_nodes=50)

# Return the validation datframwe for just the network that was used in the pipeline
df_convert_out_subset, positive_genes = myclass.alter_validation_df()

# Save a few things for checking
df_probs.to_csv(osp.join(outdir, "df_probs.tsv"), sep="\t", header=True, index=False)
df_sim_GO.to_csv(osp.join(outdir, "df_sim_GO.tsv"), sep="\t", header=True, index=False)
df_convert_out_subset.to_csv(osp.join(outdir, "df_convert_out_subset.tsv"), sep="\t", header=True, index=False)
