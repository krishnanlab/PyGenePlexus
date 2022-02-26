import os
import os.path as osp
import pathlib
import time

import numpy as np

from geneplexus import geneplexus

# The first step is for a user to load up a set of genes as a list
# this file can be found in the repo
input_genes = np.loadtxt("input_genes.txt", dtype=str, delimiter=", ")
input_genes = [item.strip("'") for item in input_genes]

# Get the data from Azure
homedir = pathlib.Path(__file__).absolute().parent
datadir = osp.join(homedir, "data")
os.makedirs(datadir, exist_ok=True)
print(f"Start downloading data and saving to: {datadir}")
for GSCs in ["GO", "DisGeNet"]:
    geneplexus.download_select_data(
        datadir,
        tasks="All",
        networks="BioGRID",
        features="Embedding",
        GSCs=GSCs,
    )
print("Done downlaoding")

# Run through the pipeline
# First initialize the geneplexus object
myclass = geneplexus.GenePlexus(datadir)

# Load the input genes into the class
myclass.load_genes(input_genes)

# Convert the input genes to Entrez
# This will return a dataframe of how the genes are in each network
df_convert_out = myclass.convert_to_Entrez()

# Set the params you want for the rest of the pipeline
myclass.set_params("BioGRID", "Embedding", "GO")

# This gets the postives and negatvies
pos_genes_in_net, negative_genes, net_genes = myclass.get_pos_and_neg_genes()

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
df_probs.to_csv("df_probs.tsv", sep="\t", header=True, index=False)
df_sim_GO.to_csv("df_sim_GO.tsv", sep="\t", header=True, index=False)
df_convert_out_subset.to_csv("df_convert_out_subset.tsv", sep="\t", header=True, index=False)
