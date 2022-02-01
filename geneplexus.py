import numpy as np
import pandas as pd
import pickle
import utls
from scipy.stats import hypergeom
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import average_precision_score
import time
from scipy.spatial.distance import cosine
import os


################################################################################################################################

# This set of functions is for running the main parts of the pipeline

'''
To Do

1. Make good repr printouts
2. Add documention of functions
3. Clear args when set_params_reintilialize is run
4. change pickle files to json (if below python 3.8 can't read the pickle files)
5. Compare all functions to what is in geneplexus_app
'''


class GenePlexus:
    
    def __init__(self,file_loc='/Users/christophermancuso/Documents/DataSets/Geneplexus_data/',network='BioGRID',
                      features='Embedding',GSC='GO'):
        self.file_loc = file_loc
        self.network = network
        self.features = features
        self.GSC = GSC
        
    def load_genes(self,input_genes):
        self.input_genes = input_genes
        
    def convert_to_Entrez(self):
        convert_IDs, df_convert_out = utls.intial_ID_convert(self.input_genes,self.file_loc)
        df_convert_out, table_summary, input_count = utls.make_validation_df(df_convert_out,self.file_loc)
        self.convert_IDs = convert_IDs
        self.df_convert_out = df_convert_out
        self.table_summary = table_summary
        self.input_count = input_count
        return self.df_convert_out
        
    def set_params(self,net_type,features,GSC):
        self.net_type = net_type
        self.features = features
        self.GSC = GSC

    def get_pos_and_neg_genes(self):
        pos_genes_in_net, genes_not_in_net, net_genes = utls.get_genes_in_network(self.file_loc,self.net_type,self.convert_IDs)
        self.pos_genes_in_net = pos_genes_in_net
        self.genes_not_in_net = genes_not_in_net
        self.net_genes = net_genes
        negative_genes = utls.get_negatives(self.file_loc,self.net_type,self.GSC,self.pos_genes_in_net)
        self.negative_genes = negative_genes
        return self.pos_genes_in_net, self.negative_genes, self.net_genes
        
    def fit_and_predict(self):
        mdl_weights, probs, avgps = utls.run_SL(self.file_loc,self.net_type,self.features,
                                                self.pos_genes_in_net,self.negative_genes,self.net_genes)
        self.mdl_weights = mdl_weights
        self.probs = probs
        self.avgps = avgps
        df_probs = utls.make_prob_df(self.file_loc,self.net_genes,self.probs,
                                     self.pos_genes_in_net,self.negative_genes)
        self.df_probs = df_probs
        return self.mdl_weights, self.df_probs, self.avgps
        
    # def make_prob_df(self):
    #     df_probs = utls.make_prob_df(self.file_loc,self.net_genes,self.probs,
    #                                  self.pos_genes_in_net,self.negative_genes)
    #     self.df_probs = df_probs
    #     return self.df_probs
        
    def make_sim_dfs(self):
        df_sim_GO, df_sim_Dis, weights_GO, weights_Dis = utls.make_sim_dfs(self.file_loc,self.mdl_weights,self.GSC,
                                                                           self.net_type,self.features)
        self.df_sim_GO = df_sim_GO
        self.df_sim_Dis = df_sim_Dis
        self.weights_GO = weights_GO
        self.weights_Dis = weights_Dis
        return self.df_sim_GO, self.df_sim_Dis, self.weights_GO, self.weights_Dis
        
    def make_small_edgelist(self,num_nodes=50):
        df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = utls.make_small_edgelist(self.file_loc,self.df_probs,
                                                                                            self.net_type,num_nodes=50)
        self.df_edge = df_edge
        self.isolated_genes = isolated_genes
        self.df_edge_sym = df_edge_sym
        self.isolated_genes_sym = isolated_genes_sym
        return self.df_edge, self.isolated_genes, self.df_edge_sym, self.isolated_genes_sym
        
    def alter_validation_df(self):
        df_convert_out_subset, positive_genes = utls.alter_validation_df(self.df_convert_out,self.table_summary,self.net_type)
        self.df_convert_out_subset = df_convert_out_subset
        self.positive_genes = positive_genes
        return self.df_convert_out_subset, self.positive_genes


################################################################################################################################
# functions specific to the webserver

# def make_graph(df_edge, df_probs):
#     max_num_genes = 1000
#     df_edge.fillna(0)
#     df_edge.columns = ['source', 'target', 'weight']
#     nodes = df_probs[0:max_num_genes]
#     nodes.rename(columns={'Entrez': 'id', 'Class-Label': 'Class'}, inplace=True)
#     nodes = nodes.astype({'id': int})
#
#     graph = {}
#     graph["nodes"] = nodes.to_dict(orient='records')
#     graph["links"] = df_edge.to_dict(orient='records')
#
#     return graph
#
# def make_template(jobname, net_type, features, GSC, avgps, df_probs, df_GO, df_dis, input_count, positive_genes, df_convert_out_subset, graph):
#     # Render the Jinja template, filling fields as appropriate
#     # return rendered HTML
#     # Find the module absolute path and locate templates
#
#     module_root = os.path.join(os.path.dirname(__file__), 'templates')
#     env = Environment(loader=FileSystemLoader(module_root))
#
#     # Find the absolute module path and the static files
#     context_menu_path = os.path.join(os.path.dirname(__file__), 'static', 'd3-v4-contextmenu.js')
#     with open(context_menu_path, 'r') as f:
#         context_menu_js = f.read()
#
#     tip_path = os.path.join(os.path.dirname(__file__), 'static', 'd3-tip.js')
#     with open(tip_path, 'r') as f:
#         d3_tip_js = f.read()
#
#     graph_path = os.path.join(os.path.dirname(__file__), 'static', 'graph.js')
#     with open(graph_path, 'r') as f:
#         graph_js = f.read()
#
#     datatable_path = os.path.join(os.path.dirname(__file__), 'static', 'datatable.js')
#     with open(datatable_path, 'r') as f:
#         datatable_js = f.read()
#
#     main_path = os.path.join(os.path.dirname(__file__), 'static', 'main.css')
#     with open(main_path, 'r') as f:
#         main_css = f.read()
#
#     graph_css_path = os.path.join(os.path.dirname(__file__), 'static', 'graph.css')
#     with open(graph_css_path, 'r') as f:
#         graph_css = f.read()
#
#     d3_tip_css_path = os.path.join(os.path.dirname(__file__), 'static', 'd3-tip.css')
#     with open(d3_tip_css_path, 'r') as f:
#         d3_tip_css = f.read()
#
#     template = env.get_template('result_base.html').render(
#         jobname=jobname,
#         network=net_type,
#         features=features,
#         negativeclass=GSC,
#         avgps=avgps,
#         input_count=input_count,
#         positive_genes=positive_genes,
#         context_menu_js=context_menu_js,
#         d3_tip_js=d3_tip_js,
#         graph_js=graph_js,
#         datatable_js=datatable_js,
#         main_css=main_css,
#         graph_css=graph_css,
#         d3_tip_css=d3_tip_css,
#         probs_table=df_probs.to_html(index=False, classes='table table-striped table-bordered" id = "probstable'),
#         go_table=df_GO.to_html(index=False,
#                                classes='table table-striped table-bordered nowrap" style="width: 100%;" id = "gotable'),
#         dis_table=df_dis.to_html(index=False, classes='table table-striped table-bordered" id = "distable'),
#         validate_results=df_convert_out_subset.to_html(index=False,
#                                               classes='table table-striped table-bordered" id = "validateresults'),
#         graph=graph)
#
#     return template
#
# def save_files(fp_save,jobname,df_probs,avgps):
#     if not os.path.exists(fp_save):
#         os.makedirs(fp_save)
#     df_probs.to_csv(fp_save+jobname+'--predictions.tsv',sep='\t',header=True,index=False)
#     np.savetxt(fp_save+jobname+'--CVvalues.txt',avgps,header='CVs (log2p)')
#
#
#
#
#
#
################################################################################################################################

# # This set of functions is for abstracting how a file is loaded
# fp_HPCC = '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/'
# def load_txtfile(file_type,file_loc,dtype_=str,net_type_=None,GSC_=None,target_set_=None):
#     if file_type == 'net_genes':
#         if file_loc == 'local':
#             output_txt = np.loadtxt('../data_backend2/Node_Orders/%s_nodelist.txt'%net_type_,dtype=dtype_)
#         elif file_loc == 'HPCC':
#             output_txt = np.loadtxt(fp_HPCC + 'data_backend2/Node_Orders/%s_nodelist.txt'%net_type_,dtype=dtype_)
#         elif file_loc == 'cloud':
#             raise ValueError('cloud is not yet implemented')
#     elif file_type == 'uni_genes':
#         if file_loc == 'local':
#             output_txt = np.loadtxt('../data_backend2/GSCs/%s_%s_universe.txt'%(GSC_,net_type_),dtype=dtype_)
#         elif file_loc == 'HPCC':
#             output_txt = np.loadtxt(fp_HPCC + 'data_backend2/GSCs/%s_%s_universe.txt'%(GSC_,net_type_),dtype=dtype_)
#         elif file_loc == 'cloud':
#             raise ValueError('cloud is not yet implemented')
#     elif file_type == 'GSC_order':
#         if file_loc == 'local':
#             output_txt = np.loadtxt('../data_backend2/CorrectionMatrices/%s_%s_Orders.txt'%(target_set_,net_type_),dtype=dtype_)
#         elif file_loc == 'HPCC':
#             output_txt = np.loadtxt(fp_HPCC + 'data_backend2/CorrectionMatrices/%s_%s_Orders.txt'%(target_set_,net_type_),dtype=dtype_)
#         elif file_loc == 'cloud':
#             raise ValueError('cloud is not yet implemented')
#     return output_txt
#
# def load_npyfile(file_type,file_loc,features_=None,net_type_=None,GSC_=None,target_set_=None):
#     if file_type == 'data':
#         if file_loc == 'local':
#             output_npy = np.load('../data_backend2/%s/%s_data.npy'%(features_,net_type_))
#         elif file_loc == 'HPCC':
#             output_npy = np.load(fp_HPCC + 'data_backend2/%s/%s_data.npy'%(features_,net_type_))
#         elif file_loc == 'cloud':
#             raise ValueError('cloud is not yet implemented')
#     elif file_type == 'cor_mat':
#         if file_loc == 'local':
#             output_npy = np.load('../data_backend2/CorrectionMatrices/%s_%s_%s_%s_CorMat.npy'%(GSC_,target_set_,net_type_,features_))
#         elif file_loc == 'HPCC':
#             output_npy = np.load(fp_HPCC + 'data_backend2/CorrectionMatrices/%s_%s_%s_%s_CorMat.npy'%(GSC_,target_set_,net_type_,features_))
#         elif file_loc == 'cloud':
#             raise ValueError('cloud is not yet implemented')
#     return output_npy
#
# def load_df(file_type,file_loc,sep_='\t',header_=None,net_type_=None):
#     if file_type == 'edgelist':
#         if file_loc == 'local':
#             if net_type_ == 'BioGRID':
#                 output_df = pd.read_csv('../data_backend2/Edgelists/%s.edg'%net_type_,sep=sep_,header=header_,names=['Node1','Node2'])
#             else:
#                 output_df = pd.read_csv('../data_backend2/Edgelists/%s.edg'%net_type_,sep=sep_,header=header_,names=['Node1','Node2','Weight'])
#         elif file_loc == 'HPCC':
#             if net_type_ == 'BioGRID':
#                 output_df = pd.read_csv(fp_HPCC + 'data_backend2/Edgelists/%s.edg'%net_type_,sep=sep_,header=header_,names=['Node1','Node2'])
#                 output_df["Weight"] = 1
#             else:
#                 output_df = pd.read_csv(fp_HPCC + 'data_backend2/Edgelists/%s.edg'%net_type_,sep=sep_,header=header_,names=['Node1','Node2','Weight'])
#         elif file_loc == 'cloud':
#             raise ValueError('cloud is not yet implemented')
#     return output_df
#
# def load_dict(file_type,file_loc,anIDtype_=None,GSC_=None,net_type_=None,target_set_=None,features_=None):
#     if file_type == 'to_Entrez':
#         if file_loc == 'local':
#             with open('../data_backend2/ID_conversion/Homo_sapiens__%s-to-Entrez__All-Mappings.pickle'%anIDtype_,'rb') as handle:
#                 output_dict = pickle.load(handle)
#         elif file_loc == 'HPCC':
#             with open(fp_HPCC + 'data_backend2/ID_conversion/Homo_sapiens__%s-to-Entrez__All-Mappings.pickle'%anIDtype_,'rb') as handle:
#                 output_dict = pickle.load(handle)
#         elif file_loc == 'cloud':
#             raise ValueError('cloud is not yet implemented')
#     elif file_type == 'good_sets':
#         if file_loc == 'local':
#             with open('../data_backend2/GSCs/%s_%s_GoodSets.pickle'%(GSC_,net_type_),'rb') as handle:
#                 output_dict = pickle.load(handle)
#         elif file_loc == 'HPCC':
#             with open(fp_HPCC + 'data_backend2/GSCs/%s_%s_GoodSets.pickle'%(GSC_,net_type_),'rb') as handle:
#                 output_dict = pickle.load(handle)
#         elif file_loc == 'cloud':
#             raise ValueError('cloud is not yet implemented')
#     elif file_type == 'Entrez_to_Symbol':
#         if file_loc == 'local':
#             with open('../data_backend2/ID_conversion/Homo_sapiens__Entrez-to-Symbol__All-Mappings.pickle','rb') as handle:
#                 output_dict = pickle.load(handle)
#         elif file_loc == 'HPCC':
#             with open(fp_HPCC + 'data_backend2/ID_conversion/Homo_sapiens__Entrez-to-Symbol__All-Mappings.pickle','rb') as handle:
#                 output_dict = pickle.load(handle)
#         elif file_loc == 'cloud':
#             raise ValueError('cloud is not yet implemented')
#     elif file_type == 'Entrez_to_Name':
#         if file_loc == 'local':
#             with open('../data_backend2/ID_conversion/Homo_sapiens__Entrez-to-Name__All-Mappings.pickle','rb') as handle:
#                 output_dict = pickle.load(handle)
#         elif file_loc == 'HPCC':
#             with open(fp_HPCC + 'data_backend2/ID_conversion/Homo_sapiens__Entrez-to-Name__All-Mappings.pickle','rb') as handle:
#                 output_dict = pickle.load(handle)
#         elif file_loc == 'cloud':
#             raise ValueError('cloud is not yet implemented')
#     elif file_type == 'weights':
#         if file_loc == 'local':
#             with open('../data_backend2/PreTrainedModels/%s_%s_%s_ModelWeights.pickle'%(target_set_,net_type_,features_),'rb') as handle:
#                 output_dict = pickle.load(handle)
#         elif file_loc == 'HPCC':
#             with open(fp_HPCC + 'data_backend2/PreTrainedModels/%s_%s_%s_ModelWeights.pickle'%(target_set_,net_type_,features_),'rb') as handle:
#                 output_dict = pickle.load(handle)
#         elif file_loc == 'cloud':
#             raise ValueError('cloud is not yet implemented')
#
#     return output_dict

