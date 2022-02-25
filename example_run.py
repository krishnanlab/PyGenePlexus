import argparse
import geneplexus
import numpy as np
import pandas as pd
import time


###################################################################################################################

# The first step is for a user to load up a set of genes as a list
# this file can be found in the repo
input_genes = np.loadtxt('input_genes.txt',dtype=str,delimiter=', ')
input_genes = [item.strip("'") for item in input_genes]

# get the data from Azure
geneplexus.download_select_data('/Users/christophermancuso/Documents/DataSets/Geneplexus_data/',
                                tasks='All',networks='BioGRID',features='Embedding',GSCs='GO')



###################################################################################################################
# run through the pipeline

myclass = geneplexus.GenePlexus('/Users/christophermancuso/Documents/DataSets/Geneplexus_data/')
# load the input genes into the class
myclass.load_genes(input_genes)
# convert the input genes to Entrez
# this will return a dataframe of how the genes are in each network
df_convert_out = myclass.convert_to_Entrez()
# Set the params you want for the rest of the pipeline
myclass.set_params('BioGRID','Embedding','GO')
# This gets the postives and negatvies
pos_genes_in_net, negative_genes, net_genes = myclass.get_pos_and_neg_genes()
# This trains the model and predcits on every gene in the network
mdl_weights, df_probs, avgps = myclass.fit_and_predict()
# The makes the tables that have the model weight similarity to other models trained on known GO and DisGeNet sets
df_sim_GO, df_sim_Dis, weights_GO, weights_Dis = myclass.make_sim_dfs()
# return an edgelist
df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = myclass.make_small_edgelist(num_nodes=50)
# return the validation datframwe for just the network that was used in the pipeline
df_convert_out_subset, positive_genes = myclass.alter_validation_df()

# save a few things for checking
df_probs.to_csv('df_probs.tsv',sep='\t',header=True,index=False)
df_sim_GO.to_csv('df_sim_GO.tsv',sep='\t',header=True,index=False)
df_convert_out_subset.to_csv('df_convert_out_subset.tsv',sep='\t',header=True,index=False)

    
###################################################################################################################
## old code for testing
    
# data = np.load('/Users/christophermancuso/Documents/DataSets/from_Azure/Embedding_BioGRID.npy')
# print(data.shape)
# print(data[0:10,0:10])
# uni_genes = np.loadtxt('/Users/christophermancuso/Documents/DataSets/from_Azure/GSC_DisGeNet_GIANT-TN_universe.txt',dtype=str)
# print(uni_genes.shape)
# print(uni_genes[0:10])