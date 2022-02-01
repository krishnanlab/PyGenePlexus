import argparse
import geneplexus
import numpy as np
import pandas as pd

# Main thing to discuss is how to handle the file_loc problem
# Also, when to save objects? Probbaly not have that be something we do for the user?

###################################################################################################################

# The first step is for a user to load up a set of genes as a list
# Need to figure out best way for doing multiple sets at once

input_genes = np.loadtxt('input_genes.txt',dtype=str,delimiter=', ')
input_genes = [item.strip("'") for item in input_genes]

###################################################################################################################

myclass = geneplexus.GenePlexus()
# load the input genes into the class
myclass.load_genes(input_genes)
# convert the input genes to Entrez
# this will return a dataframe of how the genes are in each network
df_convert_out = myclass.convert_to_Entrez()
print(df_convert_out.head())
myclass.set_params('BioGRID','Embedding','GO')
pos_genes_in_net, negative_genes, net_genes = myclass.get_pos_and_neg_genes()
mdl_weights, df_probs, avgps = myclass.fit_and_predict()
print(list(df_probs))
df_sim_GO, df_sim_Dis, weights_GO, weights_Dis = myclass.make_sim_dfs()
df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = myclass.make_small_edgelist(num_nodes=50)
df_convert_out_subset, positive_genes = myclass.alter_validation_df()
print(df_convert_out_subset.head())


    

    
