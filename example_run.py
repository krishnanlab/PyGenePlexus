import argparse
import geneplexus
import numpy as np
import pandas as pd
import time


tic = time.time()
# Main thing to discuss is how to handle the file_loc problem
# Also, when to save objects? Probbaly not have that be something we do for the user?

###################################################################################################################

# The first step is for a user to load up a set of genes as a list
# Need to figure out best way for doing multiple sets at once

input_genes = np.loadtxt('input_genes.txt',dtype=str,delimiter=', ')
input_genes = [item.strip("'") for item in input_genes]

# geneplexus.download_IDconversion_data('/Users/christophermancuso/Documents/DataSets/from_Azure/')
# geneplexus.download_all_data('/Users/christophermancuso/Documents/DataSets/from_Azure/')
geneplexus.download_select_data('/Users/christophermancuso/Documents/DataSets/from_Azure/',tasks='MachineLearning')

print('The time is took to run script is',time.time()-tic)

# data = np.load('/Users/christophermancuso/Documents/DataSets/from_Azure/Embedding_BioGRID.npy')
# print(data.shape)
# print(data[0:10,0:10])
# uni_genes = np.loadtxt('/Users/christophermancuso/Documents/DataSets/from_Azure/GSC_DisGeNet_GIANT-TN_universe.txt',dtype=str)
# print(uni_genes.shape)
# print(uni_genes[0:10])

###################################################################################################################

# myclass = geneplexus.GenePlexus()
# # load the input genes into the class
# myclass.load_genes(input_genes)
# # convert the input genes to Entrez
# # this will return a dataframe of how the genes are in each network
# df_convert_out = myclass.convert_to_Entrez()
# print(df_convert_out.head())
# myclass.set_params('BioGRID','Embedding','GO')
# pos_genes_in_net, negative_genes, net_genes = myclass.get_pos_and_neg_genes()
# mdl_weights, df_probs, avgps = myclass.fit_and_predict()
# print(list(df_probs))
# df_sim_GO, df_sim_Dis, weights_GO, weights_Dis = myclass.make_sim_dfs()
# df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = myclass.make_small_edgelist(num_nodes=50)
# df_convert_out_subset, positive_genes = myclass.alter_validation_df()
# print(df_convert_out_subset.head())


    

    
