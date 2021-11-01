import numpy as np
import pandas as pd
import pickle
import argparse
import sys
sys.path.insert(1, '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/src_backend')
import utls2
import time

parser = argparse.ArgumentParser()

parser.add_argument('-stm','--start_model',default = 0, type = int,
    help = 'The index of the target gene array to start with')
parser.add_argument('-spm','--stop_model',default = 5, type = int,
    help = 'The index of the target gene array to stop with')
parser.add_argument('-g', '--gsc', default = 'DisGeNet', type = str, 
    help = 'Which geneset collection to do')
parser.add_argument('-n','--net_type', default = 'BioGRID', type=str,
    help = 'Type of network to use')
parser.add_argument('-f','--features', default = 'Embedding', type = str, 
    help = 'Type of features to use')

args = parser.parse_args()
start_model       = args.start_model
stop_model        = args.stop_model
gsc               = args.gsc
net_type          = args.net_type
features          = args.features

tic = time.time()

with open('../data_backend2/GSCs/%s_%s_GoodSets.pickle'%(gsc,net_type),'rb') as handle:
    good_sets = pickle.load(handle)
sorted_keys = sorted(list(good_sets))

weight_dict = {}
for idx in range(start_model,stop_model,1):
    ID_tmp = sorted_keys[idx]
    weight_dict[ID_tmp] = {}
    Name_tmp = good_sets[sorted_keys[idx]]['Name']
    weight_dict[ID_tmp]['Name'] = Name_tmp
    pos_genes_tmp = good_sets[sorted_keys[idx]]['Genes']
    pos_genes_in_net, genes_not_in_net, net_genes = utls2.get_genes_in_network(pos_genes_tmp,net_type)
    negative_genes = utls2.get_negatives(pos_genes_in_net,net_type,gsc)
    mdl_weights, probs, avgps = utls2.run_SL(pos_genes_in_net,negative_genes,net_genes,net_type,features,False)
    weight_dict[ID_tmp]['Weights'] = mdl_weights
    weight_dict[ID_tmp]['PosGenes'] = pos_genes_in_net
fp_save = '/mnt/gs18/scratch/users/mancus16/GenePlexus/weight_dicts/'
with open(fp_save + '%s_%s_%s_%i_%i_weights.pickle'%(gsc,net_type,features,start_model,stop_model), 'wb') as handle:
    pickle.dump(weight_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

print()
print('The time it took the script to run is',time.time()-tic)