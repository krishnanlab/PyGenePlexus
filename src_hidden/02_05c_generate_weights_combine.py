import numpy as np
import pandas as pd
import os
import pickle
import time


tic = time.time()

# dir to put the slurm files
slurm_dir = '/mnt/gs18/scratch/users/mancus16/GenePlexus/slurms/'
fp_open = '/mnt/gs18/scratch/users/mancus16/GenePlexus/weight_dicts/'

for anet in ['BioGRID','STRING-EXP','STRING','GIANT-TN']:
    for afeat in ['Embedding','Adjacency','Influence']:
        for agsc in ['DisGeNet','GO']:
            with open('../data_backend2/GSCs/%s_%s_GoodSets.pickle'%(agsc,anet),'rb') as handle:
                good_sets = pickle.load(handle)
            if afeat == 'Embedding':
                set_size = 575
            else:
                set_size = 48
            num_bins = int(np.ceil(len(good_sets)/set_size))
            dict_final = {}
            for idx in range(num_bins):
                start_ind = idx*set_size
                stop_ind = idx*set_size + set_size
                if stop_ind > len(good_sets):
                    stop_ind = len(good_sets)
                with open(fp_open + '%s_%s_%s_%i_%i_weights.pickle'%(agsc,anet,afeat,start_ind,stop_ind),'rb') as handle:
                    dict_tmp = pickle.load(handle)
                dict_final.update(dict_tmp)
            if len(good_sets) != len(dict_final):
                print('Not size match for',anet,afeat,agsc)
            else:
                with open('../data_backend2/PreTrainedModels/%s_%s_%s_ModelWeights.pickle'%(agsc,anet,afeat), 'wb') as handle:
                    pickle.dump(dict_final, handle, protocol=pickle.HIGHEST_PROTOCOL)

