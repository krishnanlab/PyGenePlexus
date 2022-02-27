"""Generate results to be used for testing."""
import os
import os.path as osp
import pathlib
import time

import numpy as np

import geneplexus

# The first step is for a user to load up a set of genes as a list
# this file can be found in the repo
input_genes = geneplexus.util.read_gene_list("input_genes.txt")

# Set up directories
homedir = pathlib.Path(__file__).absolute().parent
datadir = osp.join(homedir, "data")
outdir = osp.join(homedir, "test", "expected_result")
os.makedirs(datadir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)

# Get the data from Azure
geneplexus.download.download_select_data(
    datadir,
    tasks="All",
    networks="BioGRID",
    features="Embedding",
    GSCs=["GO", "DisGeNet"],
)

myclass = geneplexus.GenePlexus(datadir, "BioGRID", "Embedding", "GO")
myclass.load_genes(input_genes)

myclass.convert_to_Entrez()
myclass.df_convert_out.to_csv(osp.join(outdir, "df_convert_out.tsv"), sep="\t", index=False)
# print(f"{myclass.df_convert_out=}")

myclass.get_pos_and_neg_genes()
myclass.fit_and_predict()
myclass.df_probs.to_csv(osp.join(outdir, "df_probs.tsv"), sep="\t", index=False)
# print(f"{myclass.df_probs=}")

myclass.make_sim_dfs()
myclass.df_sim_GO.to_csv(osp.join(outdir, "df_sim_GO.tsv"), sep="\t", index=False)
myclass.df_sim_Dis.to_csv(osp.join(outdir, "df_sim_Dis.tsv"), sep="\t", index=False)
# print(f"{myclass.df_sim_GO=}")
# print(f"{myclass.df_sim_Dis=}")

myclass.make_small_edgelist(num_nodes=50)
myclass.df_edge.to_csv(osp.join(outdir, "df_edge.tsv"), sep="\t", index=False)
myclass.df_edge_sym.to_csv(osp.join(outdir, "df_edge_sym.tsv"), sep="\t", index=False)
# print(f"{myclass.df_edge=}")
# print(f"{myclass.df_edge_sym=}")

myclass.alter_validation_df()
myclass.df_convert_out_subset.to_csv(
    osp.join(outdir, "df_convert_out_subset.tsv"),
    sep="\t",
    index=False,
)
# print(f"{myclass.df_convert_out_subset}")
