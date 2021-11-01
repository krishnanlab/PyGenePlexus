import numpy as np
import pandas as pd
import pickle


for net_type in ['BioGRID', 'STRING-EXP', 'STRING', 'GIANT-TN']:

    with open('../data_hidden2/embeddings/%s.emd'%net_type,'r') as f:
        for idx, line in enumerate(f):
            num_rows = int(line.strip().split()[0])
            num_cols = int(line.strip().split()[1])
            break
    data_matrix = np.zeros((num_rows,num_cols),dtype=float)
    nodelist = []
    with open('../data_hidden2/embeddings/%s.emd'%net_type,'r') as f:
        for idx, line in enumerate(f):
            if idx == 0:
                continue
            line = line.strip().split()
            nodelist.append(line[0])
            data_matrix[idx-1,:] = line[1:]
    np.savetxt('../data_backend2/Node_Orders/%s_nodelist.txt'%net_type,nodelist,fmt='%s')
    np.save('../data_backend2/Embedding/%s_data.npy'%net_type,data_matrix)
        
    