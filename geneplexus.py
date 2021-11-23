import argparse
import utls
import numpy as np


'''

1. Don't know what to do with different input types, right now hard coded in
2. How to save outputs or store them for best use in cloud
3. How to get this to run on command line (i.e. geneplexus -t validate)
4. Can we make a repo/package that is good for cloud and stand alone version

'''
### there is no logger functions in here like Doug has ####
### the way the user gene list files are read in is different ###
### the way the backend data is read in is different ###
### how to handle how webserver remembers the validate IDs then select job parameters ###
    
def read_input_file(file_path,sep=', '):
    input_genes = np.loadtxt('input_genes.txt',dtype=str,delimiter=', ')
    input_genes = [item.strip("'") for item in input_genes]
    return input_genes

def validate_genes(input_genes):
    '''
    This function is supposed to return some  info on how many genes were able to be converted
    1. This is returned on webserver, but for here maybe need something different or generate this as an output file
    '''
    convert_IDs, df_convert_out = utls.intial_ID_convert(input_genes,file_loc='HPCC')
    df_convert_out, table_summary, input_count = utls.make_validation_df(df_convert_out,file_loc='HPCC')
    print('The number of input genes is',input_count)
    print('The summary is')
    print(table_summary)
    print('The df_convert head is')
    print(df_convert_out.head())
    
def run_model(input_genes, net_type, GSC, features, jobname):

    print('0. redo get validation') ### this is in here as not sure how this is handle in webserver
    convert_IDs, df_convert_out = utls.intial_ID_convert(input_genes,file_loc='HPCC')
    df_convert_out, table_summary, input_count = utls.make_validation_df(df_convert_out,file_loc='HPCC')
    
    print('1. get_genese_in_network')
    pos_genes_in_net, genes_not_in_net, net_genes = utls.get_genes_in_network(convert_IDs,net_type,file_loc='HPCC')

    print('2. get_negatives')
    negative_genes = utls.get_negatives(pos_genes_in_net,net_type,GSC,file_loc='HPCC')

    print('3. run_SL... features=%s'%features)
    mdl_weights, probs, avgps = utls.run_SL(pos_genes_in_net,negative_genes,net_genes,
                                            net_type,features,file_loc='HPCC')

    print('4. make_prob_df...')
    df_probs, Entrez_to_Symbol = utls.make_prob_df(net_genes,probs,pos_genes_in_net,negative_genes,file_loc='HPCC')

    print('5. make_sim_dfs...')
    df_GO, df_dis, weights_dict_GO, weights_dict_Dis = utls.make_sim_dfs(mdl_weights,GSC,
                                                                         net_type,features,file_loc='HPCC')

    print('6. make_small_edgelist...')
    df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = utls.make_small_edgelist(df_probs,net_type,
                                                                                        Entrez_to_Symbol,
                                                                                        file_loc='HPCC')

    print('7. make_graph...')
    '''
    Here in webserver this max number of genes is some app.config.get variable.
    I've hard coded that varible into the function
    '''
    graph = utls.make_graph(df_edge, df_probs)

    print('8. alter_validation_df')
    df_convert_out_subset, positive_genes = utls.alter_validation_df(df_convert_out,table_summary,net_type)

    print('9. make_template...')
    '''
    Not sure how the job name is set in the webserver code
    '''
    template = utls.make_template(jobname, net_type, features, GSC, avgps, df_probs, df_GO,
                  df_dis, input_count, positive_genes, df_convert_out_subset, graph)

    
    ##########
    # need to wtite a function to save df_probs, graph, df_GO, df_dis and tempate
    #########
    

    
