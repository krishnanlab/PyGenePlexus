import numpy as np
import pandas as pd
import os
import pickle
import time
from scipy.spatial.distance import cosine



for anet in ['BioGRID','STRING-EXP','STRING','GIANT-TN']:
    with open('../data_backend/GSCs/GO_%s_GoodSets.pickle'%anet,'rb') as handle:
        GO_sets = pickle.load(handle)
    with open('../data_backend/GSCs/DisGeNet_%s_GoodSets.pickle'%anet,'rb') as handle:
        Dis_sets = pickle.load(handle)
    GO_order = list(GO_sets)
    Dis_order = list(Dis_sets)
    np.savetxt('../data_backend/CorrectionMatrices/GO_%s_Orders.txt'%anet,GO_order,fmt='%s')
    np.savetxt('../data_backend/CorrectionMatrices/DisGeNet_%s_Orders.txt'%anet,Dis_order,fmt='%s')
    for afeat in ['Embedding','Adjacency','Influence']:
        with open('../data_backend/PreTrainedModels/GO_%s_%s_ModelWeights.pickle'%(anet,afeat),'rb') as handle:
            GO_weights = pickle.load(handle)
        with open('../data_backend/PreTrainedModels/DisGeNet_%s_%s_ModelWeights.pickle'%(anet,afeat),'rb') as handle:
            Dis_weights = pickle.load(handle)
        full_info = {'GO':{'Weights':GO_weights,'Order':GO_order},'DisGeNet':{'Weights':Dis_weights,'Order':Dis_order}}
        
        tic = time.time()
        for acombo in [('GO','DisGeNet'),('DisGeNet','DisGeNet'),('GO','GO'),('DisGeNet','GO')]:
            gsc0 = acombo[0]
            weights0 = full_info[gsc0]['Weights']
            order0 = full_info[gsc0]['Order']
            gsc1 = acombo[1]
            weights1 = full_info[gsc1]['Weights']
            order1 = full_info[gsc1]['Order']
            num_rows = len(order0)
            num_cols = len(order1)
            data = np.zeros((num_rows,num_cols),dtype=float)
            for idx0 in range(num_rows):
                for idx1 in range(num_cols):
                    cos_sim = 1 - cosine(weights0[order0[idx0]]['Weights'],weights1[order1[idx1]]['Weights'])
                    data[idx0,idx1] = cos_sim
            np.save('../data_backend/CorrectionMatrices/%s_%s_%s_%s_CorMat.npy'%(gsc0,gsc1,anet,afeat),data)
            print('The time it took to do',gsc0,gsc1,anet,afeat,'is',time.time()-tic)
            
