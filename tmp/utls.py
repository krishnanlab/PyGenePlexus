import numpy as np
import pandas as pd
import pickle
from scipy.stats import hypergeom
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import average_precision_score
import time
from scipy.spatial.distance import cosine
from scipy.stats import rankdata
import os
import requests


################################################################################################################################

# This set of functions is for running the main parts of the pipeline
        
def intial_ID_convert(input_genes,file_loc):
    #load all the possible conversion dictionaries
    convert_types = ['ENSG','Symbol','ENSP','ENST']
    all_convert_dict = {}
    for anIDtype in convert_types:
        with open(file_loc+'IDconversion_Homo-sapiens_%s-to-Entrez.pickle'%anIDtype,'rb') as handle:
            convert_tmp = pickle.load(handle)
            convert_tmp = {akey.upper():convert_tmp[akey] for akey in convert_tmp}
        all_convert_dict[anIDtype] = convert_tmp

    # make some place holder arrays
    convert_IDs = [] # This will be a flat list for Entrez IDs to use as positives
    convert_out = [] # This will be a list of lists that will be used to tell user the conversions made
    for agene in input_genes:
        try:
            agene_int = int(agene)
            convert_out.append([agene_int,agene_int])
            convert_IDs.append(agene_int)
        except ValueError:
            for idx, anIDtype in enumerate(convert_types):
                if agene in all_convert_dict[anIDtype]:
                    convert_IDs = convert_IDs + all_convert_dict[anIDtype][agene]
                    convert_out.append([agene,', '.join(all_convert_dict[anIDtype][agene])])
                    break
                elif idx == len(convert_types)-1:
                    convert_out.append([agene,'Could Not be mapped to Entrez'])
    df_convert_out = pd.DataFrame(convert_out,columns=['Original_ID','ID_converted_to_Entrez'])
    df_convert_out = df_convert_out.astype({'Original_ID':str,'ID_converted_to_Entrez':str})
    # convert_IDs = convert_IDs
    # df_convert_out = df_convert_out
    return convert_IDs, df_convert_out
        
def make_validation_df(df_convert_out,file_loc):
    table_summary = []
    input_count = df_convert_out.shape[0]
    converted_genes = df_convert_out['ID_converted_to_Entrez'].to_numpy()
    for anet in ['BioGRID','STRING','STRING-EXP','GIANT-TN']:
        net_genes = np.loadtxt(file_loc+'NodeOrder_%s.txt'%anet,dtype=str)
        # net_genes = load_txtfile('net_genes',file_loc,net_type_=anet)
        df_tmp = df_convert_out[df_convert_out['ID_converted_to_Entrez'].isin(net_genes)]
        pos_genes_in_net = np.intersect1d(converted_genes,net_genes)
        table_row = {'Network': anet, 'NetworkGenes': len(net_genes), 'PositiveGenes': len(pos_genes_in_net)}
        table_summary.append(dict(table_row))
        tmp_ins = np.full(len(converted_genes),'N',dtype=str)
        tmp_ins[df_tmp.index.to_numpy()] = 'Y'
        df_convert_out['In %s?'%anet] = tmp_ins

    df_convert_out = df_convert_out.rename(columns = {'Original_ID': 'Original ID', 'ID_converted_to_Entrez': 'Entrez ID'})
    # table_summary = table_summary
    # input_count = input_count
    return df_convert_out, table_summary, input_count
    
def get_genes_in_network(file_loc,net_type,convert_IDs):
    net_genes = np.loadtxt(file_loc+'NodeOrder_%s.txt'%net_type,dtype=str)
    pos_genes_in_net = np.intersect1d(np.array(convert_IDs),net_genes)
    genes_not_in_net = np.setdiff1d(np.array(convert_IDs),net_genes)
    return pos_genes_in_net, genes_not_in_net, net_genes
    
def get_negatives(file_loc,net_type,GSC,pos_genes_in_net):
    uni_genes = np.loadtxt(file_loc+'GSC_%s_%s_universe.txt'%(GSC,net_type),dtype=str)
    with open(file_loc+'GSC_%s_%s_GoodSets.pickle'%(GSC,net_type),'rb') as handle:
        good_sets = pickle.load(handle)
    M = len(uni_genes)
    N = len(pos_genes_in_net)
    genes_to_remove = pos_genes_in_net
    for akey in good_sets:
        n = len(good_sets[akey]['Genes'])
        k = len(np.intersect1d(pos_genes_in_net,good_sets[akey]['Genes']))
        pval = hypergeom.sf(k-1, M, n, N)
        if pval < 0.05:
            genes_to_remove = np.union1d(genes_to_remove,good_sets[akey]['Genes'])
    negative_genes = np.setdiff1d(uni_genes,genes_to_remove)
    return negative_genes

def run_SL(file_loc,net_type,features,pos_genes_in_net,negative_genes,net_genes):
    pos_inds = [np.where(net_genes==agene)[0][0] for agene in pos_genes_in_net]
    neg_inds = [np.where(net_genes==agene)[0][0] for agene in negative_genes]
    data = np.load(file_loc+'Data_%s_%s.npy'%(features,net_type))
    std_scale = StandardScaler().fit(data)
    data   = std_scale.transform(data)
    Xdata = data[pos_inds+neg_inds,:]
    ydata = np.array([1]*len(pos_inds) + [0]*len(neg_inds))
    clf = LogisticRegression(max_iter=10000,solver='lbfgs',penalty='l2',C=1.0)
    clf.fit(Xdata,ydata)
    mdl_weights = np.squeeze(clf.coef_)
    probs = clf.predict_proba(data)[:,1]

    if len(pos_genes_in_net) < 15:
        avgps = [-10, -10, -10]
    else:
        avgps = []
        n_folds = 3
        skf= StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=None)
        for trn_inds, tst_inds in skf.split(Xdata,ydata):
            clf_cv = LogisticRegression(max_iter=10000,solver='lbfgs',penalty='l2',C=1.0)
            clf_cv.fit(Xdata[trn_inds],ydata[trn_inds])
            probs_cv = clf_cv.predict_proba(Xdata[tst_inds])[:,1]
            avgp = average_precision_score(ydata[tst_inds],probs_cv)
            num_tst_pos = np.sum(ydata[tst_inds])
            prior = num_tst_pos/Xdata[tst_inds].shape[0]
            log2_prior = np.log2(avgp/prior)
            avgps.append(log2_prior)
        # avgp = '{0:.2f}'.format(np.median(avgps)) # used in webserver but not for inflamation work
    return mdl_weights, probs, avgps

def make_prob_df(file_loc,net_genes,probs,pos_genes_in_net,negative_genes):
    with open(file_loc+'IDconversion_Homo-sapiens_Entrez-to-Symbol.pickle','rb') as handle:
        Entrez_to_Symbol = pickle.load(handle)
    with open(file_loc+'IDconversion_Homo-sapiens_Entrez-to-Name.pickle','rb') as handle:
        Entrez_to_Name = pickle.load(handle)
    prob_results = []
    for idx in range(len(net_genes)):
        if net_genes[idx] in pos_genes_in_net:
            class_label = 'P'
            novel_label = 'Known'
        elif net_genes[idx] in negative_genes:
            class_label = 'N'
            novel_label = 'Novel'
        else:
            class_label = 'U'
            novel_label = 'Novel'
        try:
            syms_tmp = '/'.join(Entrez_to_Symbol[net_genes[idx]]) #allows for multimapping
        except KeyError:
            syms_tmp = 'N/A'
        try:
            name_tmp = '/'.join(Entrez_to_Name[net_genes[idx]]) #allows for multimapping
        except KeyError:
            name_tmp = 'N/A'
        prob_results.append([net_genes[idx],syms_tmp,name_tmp,probs[idx],novel_label, class_label])
    df_probs = pd.DataFrame(prob_results,columns=['Entrez','Symbol','Name','Probability','Known/Novel','Class-Label'])
    df_probs = df_probs.astype({'Entrez':str,'Probability':float})
    df_probs = df_probs.sort_values(by=['Probability'],ascending=False)
    df_probs['Rank'] = rankdata(1/(df_probs['Probability'].to_numpy()+1e-9),method='min')
    return df_probs

def make_sim_dfs(file_loc,mdl_weights,GSC,net_type,features):
    dfs_out = []
    for target_set in ['GO', 'DisGeNet']:
        with open(file_loc+'PreTrainedWeights_%s_%s_%s.pickle'%(target_set,net_type,features),'rb') as handle:
            weights_dict = pickle.load(handle)
        if target_set == 'GO':
            weights_dict_GO = weights_dict
        if target_set == 'DisGeNet':
            weights_dict_Dis = weights_dict
        order = np.loadtxt(file_loc+'CorrectionMatrixOrder_%s_%s.txt'%(target_set,net_type),dtype=str)
        cor_mat = np.load(file_loc+'CorrectionMatrix_%s_%s_%s_%s.npy'%(GSC,target_set,net_type,features))
        add_row = np.zeros((1,len(order)))
        for idx, aset in enumerate(order):
            cos_sim = 1 - cosine(weights_dict[aset]['Weights'],mdl_weights)
            add_row[0,idx] = cos_sim
        cor_mat = np.concatenate((cor_mat,add_row),axis=0)
        last_row = cor_mat[-1,:]
        zq = np.maximum(0, (last_row - np.mean(last_row)) / np.std(last_row))
        zs = np.maximum(0, (last_row - np.mean(cor_mat,axis=0)) / np.std(cor_mat,axis=0))
        z = np.sqrt(zq**2 + zs**2)
        results_tmp = []
        for idx2, termID_tmp in enumerate(order):
            ID_tmp = termID_tmp
            Name_tmp = weights_dict[termID_tmp]['Name']
            z_tmp = z[idx2]
            results_tmp.append([ID_tmp,Name_tmp,z_tmp])
        df_tmp = pd.DataFrame(results_tmp,columns=['ID','Name','Similarity']).sort_values(by=['Similarity'],ascending=False)
        df_tmp['Rank'] = rankdata(1/(df_tmp['Similarity'].to_numpy()+1e-9),method='min')
        dfs_out.append(df_tmp)
    return dfs_out[0], dfs_out[1], weights_dict_GO, weights_dict_Dis

def make_small_edgelist(file_loc,df_probs,net_type,num_nodes=50):
    # This will set the max number of genes to look at to a given number
    if net_type == 'BioGRID':
        df_edge = pd.read_csv(file_loc+'Edgelist_%s.edg'%net_type,sep='\t',header=None,names=['Node1','Node2'])
        df_edge['Weight'] = 1
    else:
        df_edge = pd.read_csv(file_loc+'Edgelist_%s.edg'%net_type,sep='\t',header=None,names=['Node1','Node2','Weight'])
    df_edge = df_edge.astype({'Node1':str,'Node2':str})
    top_genes = df_probs['Entrez'].to_numpy()[0:num_nodes]
    df_edge = df_edge[(df_edge['Node1'].isin(top_genes)) & (df_edge['Node2'].isin(top_genes))]
    genes_in_edge = np.union1d(df_edge['Node1'].unique(),df_edge['Node2'].unique())
    isolated_genes = np.setdiff1d(top_genes,genes_in_edge)
    with open(file_loc+'IDconversion_Homo-sapiens_Entrez-to-Symbol.pickle','rb') as handle:
        Entrez_to_Symbol = pickle.load(handle)
    replace_dict = {}
    for agene in genes_in_edge:
        try:
            syms_tmp = '/'.join(Entrez_to_Symbol[agene]) #allows for multimapping
        except KeyError:
            syms_tmp = 'N/A'
        replace_dict[agene] = syms_tmp
    df_edge_sym = df_edge.replace(to_replace=replace_dict)
    # make smae network as above just with gene symbols instead of entrez IDs
    isolated_genes_sym = []
    for agene in isolated_genes:
        try:
            syms_tmp = '/'.join(Entrez_to_Symbol[agene]) #allows for multimapping
        except KeyError:
            syms_tmp = 'N/A'
        isolated_genes_sym.append(syms_tmp)
    return df_edge, isolated_genes, df_edge_sym, isolated_genes_sym

def alter_validation_df(df_convert_out,table_summary,net_type):
    df_convert_out_subset = df_convert_out[['Original ID','Entrez ID','In %s?'%net_type]]
    network = next((item for item in table_summary if item['Network'] == net_type), None)
    positive_genes = network.get("PositiveGenes")
    return df_convert_out_subset, positive_genes

################################################################################################################################

# functions for downloading data

def download_from_azure(fp_data,files_to_do):
    for afile in files_to_do:
        if os.path.exists(fp_data+afile):
            print('The following file already exsists so skipping download: ', fp_data+afile)
            print() 
        else:   
            FN_Azure = 'https://mancusogeneplexusstorage.blob.core.windows.net/mancusoblob2/%s'%afile
            print('Downloading the follwing file: ', FN_Azure)
            print()
            r = requests.get(FN_Azure)
            open(fp_data+afile, 'wb').write(r.content)
    

def make_download_options_lists(tasks,networks,features,GSCs):
    if isinstance(tasks, str):
        if tasks == 'All':
            tasks = ['IDconversion','MachineLearning','Similarities','NetworkGraph']
        else:
            tasks = [tasks]
    if isinstance(networks, str):
        if networks == 'All':
            networks = ['BioGRID','STRING','STRING-EXP','GIANT-TN']
        else:
            networks = [networks]
    if isinstance(features, str):
        if features == 'All':
            features = ['Adjacency','Embedding','Influence']
        else:
            features = [features]
    if isinstance(GSCs, str):
        if GSCs == 'All':
            GSCs = ['GO','DisGeNet']
        else:
            GSCs = [GSCs]
    return tasks, networks, features, GSCs
    
def get_IDconversion_filenames():
    files_to_do = []
    with open('data_filenames.txt','r') as f:
        for line in f:
            line = line.strip()
            if ('IDconversion' in line) or ('NodeOrder' in line):
                files_to_do.append(line)
    return files_to_do
    
def get_MachineLearning_filenames(networks,GSCs,features):
    files_to_do = []
    with open('data_filenames.txt','r') as f:
        for line in f:
            line = line.strip()
            if 'NodeOrder' in line:
                net_tmp = line.split('Order_')[-1].split('.tx')[0]
                if net_tmp in networks:
                    files_to_do.append(line)
            if ('universe.txt' in line) or ('GoodSets.pickle' in line):
                net_tmp = line.split('_')[2]
                GSC_tmp = line.split('_')[1]
                if (net_tmp in networks) and (GSC_tmp in GSCs):
                    files_to_do.append(line)
            if 'Data_' in line:
                feature_tmp = line.split('_')[1]
                net_tmp = line.split('_')[2].split('.n')[0]
                if (feature_tmp in features) and (net_tmp in networks):
                    files_to_do.append(line)
            if ('Entrez-to-Name' in line) or ('Entrez-to-Symbol' in line):
                files_to_do.append(line)
    return files_to_do
    
def get_Similarities_filenames(networks, features, GSCs):
    files_to_do = []
    with open('data_filenames.txt','r') as f:
        for line in f:
            line = line.strip()
            if 'CorrectionMatrix_' in line:
                feature_tmp = line.split('_')[-1].split('.n')[0]
                net_tmp = line.split('_')[3]
                GSC_tmp = line.split('_')[1]
                if (net_tmp in networks) and (feature_tmp in features) and (GSC_tmp in GSCs):
                    files_to_do.append(line)
            if 'CorrectionMatrixOrder' in line:
                GSC_tmp = line.split('_')[1]
                net_tmp = line.split('_')[2].split('.t')[0]
                if (net_tmp in networks) and (GSC_tmp in GSCs):
                    files_to_do.append(line)
            if 'PreTrainedWeights' in line:
                net_tmp = line.split('_')[2]
                feature_tmp = line.split('_')[-1].split('.pic')[0]
                if (net_tmp in networks) and (feature_tmp in features):
                    files_to_do.append(line)
    return files_to_do
    
def get_NetworkGraph_filenames(networks):
    files_to_do = ['IDconversion_Homo-sapiens_Entrez-to-Symbol.pickle']
    with open('data_filenames.txt','r') as f:
        for line in f:
            line = line.strip()
            if 'Edgelist' in line:
                net_tmp = line.split('_')[-1].split('.ed')[0]
                if net_tmp in networks:
                    files_to_do.append(line)
    return files_to_do
    

################################################################################################################################

# functions just for webserver

# def make_graph(df_edge, df_probs):
#     num_nodes = 1000
#     df_edge.fillna(0)
#     df_edge.columns = ['source', 'target', 'weight']
#     nodes = df_probs[0:num_nodes]
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

# This set of functions is for abstracting how a file is loaded
fp_HPCC = '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/'
def load_txtfile(file_type,file_loc,dtype_=str,net_type_=None,GSC_=None,target_set_=None):
    if file_type == 'net_genes':
        output_txt = np.loadtxt(file_loc+'%s_nodelist.txt'%net_type_,dtype=dtype_)
    elif file_type == 'uni_genes':
        if file_loc == 'local':
            output_txt = np.loadtxt('../data_backend2/GSCs/%s_%s_universe.txt'%(GSC_,net_type_),dtype=dtype_)
        elif file_loc == 'HPCC':
            output_txt = np.loadtxt(fp_HPCC + 'data_backend2/GSCs/%s_%s_universe.txt'%(GSC_,net_type_),dtype=dtype_)
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')
    elif file_type == 'GSC_order':
        if file_loc == 'local':
            output_txt = np.loadtxt('../data_backend2/CorrectionMatrices/%s_%s_Orders.txt'%(target_set_,net_type_),dtype=dtype_)
        elif file_loc == 'HPCC':
            output_txt = np.loadtxt(fp_HPCC + 'data_backend2/CorrectionMatrices/%s_%s_Orders.txt'%(target_set_,net_type_),dtype=dtype_)
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')
    return output_txt

def load_npyfile(file_type,file_loc,features_=None,net_type_=None,GSC_=None,target_set_=None):
    if file_type == 'data':
        if file_loc == 'local':
            output_npy = np.load('../data_backend2/%s/%s_data.npy'%(features_,net_type_))
        elif file_loc == 'HPCC':
            output_npy = np.load(fp_HPCC + 'data_backend2/%s/%s_data.npy'%(features_,net_type_))
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')
    elif file_type == 'cor_mat':
        if file_loc == 'local':
            output_npy = np.load('../data_backend2/CorrectionMatrices/%s_%s_%s_%s_CorMat.npy'%(GSC_,target_set_,net_type_,features_))
        elif file_loc == 'HPCC':
            output_npy = np.load(fp_HPCC + 'data_backend2/CorrectionMatrices/%s_%s_%s_%s_CorMat.npy'%(GSC_,target_set_,net_type_,features_))
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')
    return output_npy

def load_df(file_type,file_loc,sep_='\t',header_=None,net_type_=None):
    if file_type == 'edgelist':
        if file_loc == 'local':
            if net_type_ == 'BioGRID':
                output_df = pd.read_csv('../data_backend2/Edgelists/%s.edg'%net_type_,sep=sep_,header=header_,names=['Node1','Node2'])
            else:
                output_df = pd.read_csv('../data_backend2/Edgelists/%s.edg'%net_type_,sep=sep_,header=header_,names=['Node1','Node2','Weight'])
        elif file_loc == 'HPCC':
            if net_type_ == 'BioGRID':
                output_df = pd.read_csv(fp_HPCC + 'data_backend2/Edgelists/%s.edg'%net_type_,sep=sep_,header=header_,names=['Node1','Node2'])
                output_df["Weight"] = 1
            else:
                output_df = pd.read_csv(fp_HPCC + 'data_backend2/Edgelists/%s.edg'%net_type_,sep=sep_,header=header_,names=['Node1','Node2','Weight'])
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')
    return output_df

def load_dict(file_type,file_loc,anIDtype_=None,GSC_=None,net_type_=None,target_set_=None,features_=None):
    if file_type == 'to_Entrez':
        with open(file_loc+'Homo_sapiens__%s-to-Entrez__All-Mappings.pickle'%anIDtype_,'rb') as handle:
            output_dict = pickle.load(handle)
    elif file_type == 'good_sets':
        if file_loc == 'local':
            with open('../data_backend2/GSCs/%s_%s_GoodSets.pickle'%(GSC_,net_type_),'rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'HPCC':
            with open(fp_HPCC + 'data_backend2/GSCs/%s_%s_GoodSets.pickle'%(GSC_,net_type_),'rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')
    elif file_type == 'Entrez_to_Symbol':
        if file_loc == 'local':
            with open('../data_backend2/ID_conversion/Homo_sapiens__Entrez-to-Symbol__All-Mappings.pickle','rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'HPCC':
            with open(fp_HPCC + 'data_backend2/ID_conversion/Homo_sapiens__Entrez-to-Symbol__All-Mappings.pickle','rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')
    elif file_type == 'Entrez_to_Name':
        if file_loc == 'local':
            with open('../data_backend2/ID_conversion/Homo_sapiens__Entrez-to-Name__All-Mappings.pickle','rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'HPCC':
            with open(fp_HPCC + 'data_backend2/ID_conversion/Homo_sapiens__Entrez-to-Name__All-Mappings.pickle','rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')
    elif file_type == 'weights':
        if file_loc == 'local':
            with open('../data_backend2/PreTrainedModels/%s_%s_%s_ModelWeights.pickle'%(target_set_,net_type_,features_),'rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'HPCC':
            with open(fp_HPCC + 'data_backend2/PreTrainedModels/%s_%s_%s_ModelWeights.pickle'%(target_set_,net_type_,features_),'rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')

    return output_dict

