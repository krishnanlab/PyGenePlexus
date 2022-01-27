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

# # in the webserver a user can first validate their genes and this section captures that
# # uses input genes from above
#
# convert_IDs, df_convert_out = geneplexus.intial_ID_convert(input_genes,file_loc='HPCC')
# df_convert_out, table_summary, input_count = geneplexus.make_validation_df(df_convert_out,file_loc='HPCC')
# print('The number of input genes is',input_count)
# print('The summary is')
# print(table_summary)
# print('The df_convert head is')
# print(df_convert_out.head())

###################################################################################################################

# This steps through all the functions used in the pipeline
# In addtion to below choices it uses input_genes from above

net_type = 'BioGRID'
GSC = 'GO'
features = 'Embedding'
jobname = 'mytest1'
fp_save = 'results/'


print('1. Making validation datframes')
convert_IDs, df_convert_out = geneplexus.intial_ID_convert(input_genes,file_loc='HPCC')
df_convert_out, table_summary, input_count = geneplexus.make_validation_df(df_convert_out,file_loc='HPCC')

print('2. Finding genes in network and which are positives')
pos_genes_in_net, genes_not_in_net, net_genes = geneplexus.get_genes_in_network(convert_IDs,net_type,file_loc='HPCC')

print('3. Finding which genes should be negatives')
negative_genes = geneplexus.get_negatives(pos_genes_in_net,net_type,GSC,file_loc='HPCC')

print('4. Doing the machine learning part')
mdl_weights, probs, avgps = geneplexus.run_SL(pos_genes_in_net,negative_genes,net_genes,
                                              net_type,features,file_loc='HPCC')

print('5. Making dategrame of the predictions')
df_probs, Entrez_to_Symbol = geneplexus.make_prob_df(net_genes,probs,pos_genes_in_net,negative_genes,file_loc='HPCC')


print('6. Making dataframes for the similarity tables')
df_GO, df_dis, weights_dict_GO, weights_dict_Dis = geneplexus.make_sim_dfs(mdl_weights,GSC,
                                                                     net_type,features,file_loc='HPCC')

print('7. Making edgleist for the graph')
df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = geneplexus.make_small_edgelist(df_probs,net_type,
                                                                                          Entrez_to_Symbol,
                                                                                          file_loc='HPCC')

print('8. Making a D3 readable graph object')
'''
Here in webserver this max number of genes is some app.config.get variable.
I've hard coded that varible into the function
'''
graph = geneplexus.make_graph(df_edge, df_probs)

print('9. Making output validation dataframe')
# This differs from step one in that it is only for the selected network
df_convert_out_subset, positive_genes = geneplexus.alter_validation_df(df_convert_out,table_summary,net_type)

print('10. Making html template')
'''
Not sure how the job name is set in the webserver code
'''
template = geneplexus.make_template(jobname, net_type, features, GSC, avgps, df_probs, df_GO,
                                    df_dis, input_count, positive_genes, df_convert_out_subset, graph)

print('11. Saving some things')
# right now only saving a few objects
geneplexus.save_files(fp_save,jobname,df_probs,avgps)


    

    
