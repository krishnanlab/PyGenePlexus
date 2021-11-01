import numpy as np
import pandas as pd
import pickle
import argparse
import sys
sys.path.insert(1, '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/src_backend')
import utls
import time

'''
This will add post genes to the weight_dicts but should be done in orginal script next time
'''

for anet in ['BioGRID','STRING','STRING-EXP','GIANT-TN']:
    for agsc in ['GO','DisGeNet']:
        for afeat in ['Embedding','Adjacency','Influence']:
            with open('../data_backend/GSCs/%s_%s_GoodSets.pickle'%(agsc,anet),'rb') as handle:
                good_sets = pickle.load(handle)
            sorted_keys = sorted(list(good_sets))
            with open('../data_backend/PreTrainedModels/%s_%s_%s_ModelWeights.pickle'%(agsc,anet,afeat),'rb') as handle:
                weight_dict = pickle.load(handle)
            for akey in weight_dict:
                weight_dict[akey]['PosGenes'] = good_sets[akey]['Genes']
            with open('../data_backend/PreTrainedModels/%s_%s_%s_ModelWeights.pickle'%(agsc,anet,afeat),'wb') as handle:
                pickle.dump(weight_dict,handle,protocol=pickle.HIGHEST_PROTOCOL)
            
