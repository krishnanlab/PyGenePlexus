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


################################################################################################################################

# This set of functions is for running the main parts of the pipeline

def intial_ID_convert(input_genes,file_loc='local'):
    #load all the possible conversion dictionaries 
    convert_types = ['ENSG','Symbol','ENSP','ENST']
    all_convert_dict = {}
    for anIDtype in convert_types:
        convert_tmp = load_dict('to_Entrez',file_loc,anIDtype_=anIDtype)
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
    return convert_IDs, df_convert_out
    
def make_validation_df(df_convert_out,file_loc='local'):
    converted_genes = df_convert_out['ID_converted_to_Entrez'].tolist()
    summary_results = [] 
    for anet in ['BioGRID','STRING','STRING-EXP','GIANT-TN']:
        net_genes = load_txtfile('net_genes',file_loc,net_type_=anet)
        in_net = []
        in_net_count = 0
        for agene in converted_genes:
            if agene in net_genes:
                in_net.append('Y')
                in_net_count = in_net_count + 1
            else:
                in_net.append('N')
        summary_results.append([anet,len(net_genes),in_net_count])
        df_convert_out['In %s?'%anet] = in_net
    df_summary = pd.DataFrame(summary_results,columns=['Network','Num. of Network Genes', 'Num. Positive Genes']) 
    return df_convert_out, df_summary
        
def get_genes_in_network(convert_IDs,net_type,file_loc='local'):
    net_genes = load_txtfile('net_genes',file_loc,net_type_=net_type)
    pos_genes_in_net = np.intersect1d(np.array(convert_IDs),net_genes)
    genes_not_in_net = np.setdiff1d(np.array(convert_IDs),net_genes)
    return pos_genes_in_net, genes_not_in_net, net_genes
    
def get_negatives(pos_genes_in_net,net_type,GSC,file_loc='local'):
    uni_genes = load_txtfile('uni_genes',file_loc,net_type_=net_type,GSC_=GSC)
    good_sets = load_dict('good_sets',file_loc,GSC_=GSC,net_type_=net_type)
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
    
def run_SL(pos_genes_in_net,negative_genes,net_genes,net_type,features,CV,file_loc='local'):
    pos_inds = [np.where(net_genes==agene)[0][0] for agene in pos_genes_in_net]
    neg_inds = [np.where(net_genes==agene)[0][0] for agene in negative_genes]
    data = load_npyfile('data',file_loc,features_=features,net_type_=net_type)
    
    std_scale = StandardScaler().fit(data)
    data   = std_scale.transform(data)
    Xdata = data[pos_inds+neg_inds,:]
    ydata = np.array([1]*len(pos_inds) + [0]*len(neg_inds))
    clf = LogisticRegression(max_iter=10000,solver='lbfgs',penalty='l2',C=1.0)
    clf.fit(Xdata,ydata)
    mdl_weights = np.squeeze(clf.coef_)
    probs = clf.predict_proba(data)[:,1]
    
    avgps = []
    n_folds = 5
    if CV == 'doCV': # add something here if number of positives is less than the folds
        if len(pos_genes_in_net) < n_folds:
            pass
        else:
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
    return mdl_weights, probs, avgps
    
def make_prob_df(net_genes,probs,pos_genes_in_net,negative_genes,file_loc='local'):
    Entrez_to_Symbol = load_dict('Entrez_to_Symbol',file_loc)
    Entrez_to_Name = load_dict('Entrez_to_Name',file_loc)
    prob_results = []
    for idx in range(len(net_genes)):
        if net_genes[idx] in pos_genes_in_net:
            class_label = 'P'
        elif net_genes[idx] in negative_genes:
            class_label = 'N'
        else:
            class_label = 'U'
        try:
            syms_tmp = '/'.join(Entrez_to_Symbol[net_genes[idx]]) #allows for multimapping
        except KeyError:
            syms_tmp = 'N/A'
        try:
            name_tmp = '/'.join(Entrez_to_Name[net_genes[idx]]) #allows for multimapping
        except KeyError:
            name_tmp = 'N/A'
        prob_results.append([net_genes[idx],syms_tmp,name_tmp,probs[idx],class_label])
    df_probs = pd.DataFrame(prob_results,columns=['Entrez','Symbol','Name','Probability','Class-Label'])
    df_probs = df_probs.astype({'Entrez':str,'Probability':float})
    df_probs = df_probs.sort_values(by=['Probability'],ascending=False)
    return df_probs, Entrez_to_Symbol
    
def make_sim_dfs(mdl_weights,GSC,net_type,features,file_loc='local'):
    dfs_out = []
    for target_set in ['GO', 'DisGeNet']:
        weights_dict = load_dict('weights',file_loc,net_type_=net_type,target_set_=target_set,features_=features)
        if target_set == 'GO':
            weights_dict_GO = weights_dict
        if target_set == 'DisGeNet':
            weights_dict_Dis = weights_dict
        order = load_txtfile('GSC_order',file_loc,net_type_=net_type,target_set_=target_set)
        cor_mat = load_npyfile('cor_mat',file_loc,GSC_=GSC,target_set_=target_set,net_type_=net_type,features_=features)
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
        dfs_out.append(df_tmp)
    return dfs_out[0], dfs_out[1], weights_dict_GO, weights_dict_Dis
        
    
    
def make_small_edgelist(df_probs,net_type,Entrez_to_Symbol,file_loc='local'):
    # This will set the max number of genes to look at to a given number
    max_num_genes = 1000
    df_edge = load_df('edgelist',file_loc,net_type_=net_type)
    df_edge = df_edge.astype({'Node1':str,'Node2':str})
    top_genes = df_probs['Entrez'].to_numpy()[0:max_num_genes]
    df_edge = df_edge[(df_edge['Node1'].isin(top_genes)) & (df_edge['Node2'].isin(top_genes))]
    genes_in_edge = np.union1d(df_edge['Node1'].unique(),df_edge['Node2'].unique())
    isolated_genes = np.setdiff1d(top_genes,genes_in_edge)
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
    
    
################################################################################################################################

# This set of functions is for abstracting how a file is loaded
fp_HPCC = '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/'
def load_txtfile(file_type,file_loc,dtype_=str,net_type_=None,GSC_=None,target_set_=None):
    if file_type == 'net_genes':
        if file_loc == 'local':
            output_txt = np.loadtxt('../data_backend2/Node_Orders/%s_nodelist.txt'%net_type_,dtype=dtype_)
        elif file_loc == 'HPCC':
            output_txt = np.loadtxt(fp_HPCC + 'data_backend2/Node_Orders/%s_nodelist.txt'%net_type_,dtype=dtype_)
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')
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
            else:
                output_df = pd.read_csv(fp_HPCC + 'data_backend2/Edgelists/%s.edg'%net_type_,sep=sep_,header=header_,names=['Node1','Node2','Weight'])
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')
    return output_df
    
def load_dict(file_type,file_loc,anIDtype_=None,GSC_=None,net_type_=None,target_set_=None,features_=None):
    if file_type == 'to_Entrez':
        if file_loc == 'local':
            with open('../data_backend2/ID_conversion/Homo_sapiens__%s-to-Entrez__All-Mappings.pickle'%anIDtype_,'rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'HPCC':
            with open(fp_HPCC + 'data_backend2/ID_conversion/Homo_sapiens__%s-to-Entrez__All-Mappings.pickle'%anIDtype_,'rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'cloud':
            raise ValueError('cloud is not yet implemented')
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
            



# load convert Entrez to Symbol (dict)
# load weights (dict)


if __name__ == '__main__':

    #################### This section has some inputs ##############################################
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--net_type',
                        default = 'BioGRID',
                        type = str,
                        help = 'options are BioGRID, STRING-EXP, STRING, GIANT-TN')
    parser.add_argument('-f','--features',
                        default = 'Embedding',
                        type = str,
                        help = 'options are Embedding, Adjacency, Influence')
    parser.add_argument('-g','--GSC',
                        default = 'GO',
                        type = str,
                        help = 'options are GO, DisGeNet')
    parser.add_argument('-t','--term_IDnum',
                        default = 2000,
                        type = int,
                        help = 'ID to do')
    parser.add_argument('-i','--term_IDname',
                        default = 'GO:0030217',
                        type = str,
                        help = 'Id Name to do')

    args = parser.parse_args()
    net_type = args.net_type
    features = args.features
    term_IDnum = args.term_IDnum
    term_IDname = args.term_IDname
    GSC = args.GSC
    
    CV = False # options are False or True, wether to do k_fold cross validation

    print() # this print is to make new line for weird slurm output reading
    
    # this is a list of Entrez geens for some known set
    input_genes_Entrez = ['6457', '7037', '57403', '3134', '50807', '93343', '11311', '8766', '5584', '137492', '998', 
                   '30011', '5337', '3312', '155', '10015', '55738', '57132', '153', '116986', '163', '11267', 
                   '1950', '3559', '6714', '84249', '2066', '29924', '1213', '30846', '84612', '440073', '2060', 
                   '3303', '3561', '9101', '51160', '56904', '3304', '23527', '5878', '3560', '7189', '3949', 
                   '92421', '26286', '5979', '9922', '11031', '116983', '2261', '9230', '5867', '64145', 
                   '867', '57154', '84313', '3577', '116987', '10617', '1436', '200576', '83737', '23396', '3310', '5590', '3133', '382', '6456', 
                   '30845', '868', '2264', '5868', '84440', '116984', '5869', '23624', '22841', '161', 
                   '23096', '5338', '652614', '84552', '51028', '55616', '9829', '3815', '29082', '9135', '23362', '9146', '128866', '156', 
                   '8218', '89853', '154', '64744', '9525', '84364', '9727', '23550', '8853', '1956', '8395', '6455', '64411', 
                   '5156', '51100', '8027', '408', '3305', '51534', '2868', '9744', '3106', '51652', '3265', '27243', '10938', 
                   '60682', '157', '26056', '10059', '2321', '80230', '1173', '1175', '160', '3306', '3135', '1234', '2149', 
                   '8411', '3791', '51510', '23327', '409', '11059', '3579', '27183', '8396', '1601', '1211', '3480', 
                   '9815', '26119', '64750', '26052', '4914', '25978', '8394', '1212', '30844', '131890', '79720', 
                   '7251', '50855', '116985', '5662', '2870', '10193', '1785', '155382', '652799', '22905', '3105', 
                   '55048', '10254', '55040', '7852', '1759', '4193', '2869', '2065', '6011', '4734', '28964', 
                   '4233', '80223', '79643', '3107', '2263', '56288']
    # Make a new list that has the IDs converted to ENSG
    with open('../data_backend2/ID_conversion/Homo_sapiens__Entrez-to-ENSG__All-Mappings.pickle','rb') as handle:
        convert_tmp = pickle.load(handle)
    input_genes_ENSG = []
    for agene in input_genes_Entrez:
        if agene in convert_tmp:
            input_genes_ENSG = input_genes_ENSG + convert_tmp[agene]
    # this section will read the genes in from the text file
    fp_txt = '/mnt/research/compbio/krishnanlab/tmp/to_Chris_from_Chris/input_genes.txt'
    input_genes_txtfile = np.loadtxt(fp_txt,dtype=str,delimiter=', ')
    input_genes_txtfile = [item.strip("'") for item in input_genes_txtfile]
    # # section picks a geneset by termID_num
    # with open('../data_backend2/PreTrainedModels/%s_%s_%s_ModelWeights.pickle'%(GSC,net_type,features),'rb') as handle:
    #     term_dict_ = pickle.load(handle)
    # term_dict_keys = sorted(list(term_dict_))
    # print('The name of the set is',term_dict_[term_dict_keys[term_IDnum]]['Name'], 'and the list index is',term_IDnum)
    # input_genes_term_IDnum = term_dict_[term_dict_keys[term_IDnum]]['PosGenes']
    # section picks a geneset by term_IDname
    with open('../data_backend2/PreTrainedModels/%s_%s_%s_ModelWeights.pickle'%(GSC,net_type,features),'rb') as handle:
        term_dict_ = pickle.load(handle)
    print('The name of the set is',term_dict_[term_IDname]['Name'], 'and the ID is',term_IDname)
    input_genes_term_IDname = term_dict_[term_IDname]['PosGenes']
    # this sets which gene type to use
    input_genes = input_genes_term_IDname # options now are input_genes_ENSG or input_genes_Entrez but for real anythong the user uploads
     

        
    # If using Adjacency or Influence the machine learning takes somewhere around 40 seconds, Embedding it is like 2 seconds
    # If net_type is GIANT-TN the step to make a smaller edgelist takes like 2 minutes
    
    ######## some old code ##############
    # net_type = 'BioGRID' # options are 'BioGRID','STRING-EXP','STRING','GIANT-TN' (GIANT-TN is the big one, BioGRID is a small one)
    # GSC = 'GO' #options are 'GO', 'DisGeNet
    # features = 'Embedding' # options are 'Embedding','Adjacency','Influence'
    # term_IDnum = 2000 # This picks out a nknow set for testing the pipeline
    # this is used to get term names fo benchmarking
    # for GSC in ['GO','DisGeNet']:
    #     cnt = 0
    #     for net_type in ['BioGRID','STRING','STRING-EXP','GIANT-TN']:
    #         for features in ['Embedding','Adjacency','Influence']:
    #             # this section will make a set based on gmt
    #             with open('../data_backend2/PreTrainedModels/%s_%s_%s_ModelWeights.pickle'%(GSC,net_type,features),'rb') as handle:
    #                 term_dict_ = pickle.load(handle)
    #             # with open('../data_backend2/GSCs/%s_%s_GoodSets.pickle'%(GSC,net_type),'rb') as handle:
    #             #     term_dict_ = pickle.load(handle)
    #             if cnt == 0:
    #                 all_terms = np.array(list(term_dict_))
    #                 cnt = 1
    #             else:
    #                 all_terms = np.intersect1d(all_terms,np.array(list(term_dict_)))
    #     print(len(all_terms))
    #     print(np.random.choice(all_terms,3))
            
    ################## This section runs the pipeline ###############################
    # for net_type in ['BioGRID','STRING','STRING-EXP','GIANT-TN']:
    #     for features in ['Embedding','Adjacency','Influence']:
    #         for GSC in ['GO','DisGeNet']:
    print(net_type,features,GSC)
    tic = time.time()
    convert_IDs, df_convert_out = intial_ID_convert(input_genes) # converts user genes to Entrez, df_convert_out could be a downloadable file
    pos_genes_in_net, genes_not_in_net, net_genes = get_genes_in_network(convert_IDs,net_type) #genes_not_in_net could be a donloadable file
    df_convert_out = make_validation_df(df_convert_out,pos_genes_in_net)
    # print(df_convert_out.head())
    negative_genes = get_negatives(pos_genes_in_net,net_type,GSC)
    print('The number of seconds it took to do up to the SL part is',time.time()-tic)
    tic = time.time()
    mdl_weights, probs, avgps = run_SL(pos_genes_in_net,negative_genes,net_genes,net_type,features,CV)
    print('The number of seconds it took to do the SL part is',time.time()-tic)
    tic = time.time()
    df_probs, Entrez_to_Symbol = make_prob_df(net_genes,probs,pos_genes_in_net,negative_genes) # df_ptobs is file that will be displayed on webserver
    df_GO, df_dis, weights_dict_GO, weights_dict_Dis = make_sim_dfs(mdl_weights,GSC,net_type,features) # both of these dfs will be displaed on the webserver
    print('The number of seconds it took to generate the main output dfs is',time.time()-tic)
    tic = time.time()
    # The point of below function is too take the full edgelist and make much smaller by using on top N genes from prediction step
    # This takes like 2 minutes on GIANT-TN so it is mostly here so if user changes number of nodes later it should be quick using
    # this parsed down edgelist
    df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = make_small_edgelist(df_probs,net_type,Entrez_to_Symbol)
    print('The number of seconds it took to make the small network is',time.time()-tic)
    print()
    # print(df_probs.head())
    # print(df_GO.head())
    # print(df_dis.head())

    # print(list(weights_dict_GO['GO:0000002']))
    # print(list(weights_dict_Dis['DOID:0001816']))
    # fp_save = '/mnt/research/compbio/krishnanlab/tmp/to_Chris_from_Chris/'
    # df_probs.to_csv(fp_save + 'gene_probs.tsv',sep='\t',header=True,index=False)

    
    
    
