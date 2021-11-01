import numpy as np
import pandas as pd
import pickle
import time
import scipy.linalg as sla


for net_type in ['BioGRID', 'STRING-EXP', 'STRING', 'GIANT-TN']:
    print(net_type)
    tic = time.time()
    nodelist = np.loadtxt('../data_backend2/Node_Orders/%s_nodelist.txt'%net_type,dtype=str)
    # make node to ind dict
    node_to_ind = {}
    for idx, anode in enumerate(nodelist):
        node_to_ind[anode] = idx
    # make an exmpty array
    adj_mat = np.zeros((len(nodelist),len(nodelist)),dtype=float)
    with open('../data_backend2/Edgelists/%s.edg'%net_type,'r') as f:
        for idx, line in enumerate(f):
            line = line.strip().split('\t')
            if (line[0] not in node_to_ind) or (line[1] not in node_to_ind):
                raise KeyError('Nodes not matching')
            if net_type == 'BioGRID':
                adj_mat[node_to_ind[line[0]],node_to_ind[line[1]]] = 1.
                adj_mat[node_to_ind[line[1]],node_to_ind[line[0]]] = 1.
            else:
                adj_mat[node_to_ind[line[0]],node_to_ind[line[1]]] = float(line[2])
                adj_mat[node_to_ind[line[1]],node_to_ind[line[0]]] = float(line[2])
    print('Time it took to make adj_mat',time.time()-tic)
    np.save('../data_backend2/Adjacency/%s_data.npy'%net_type,adj_mat)
    tic = time.time()
    adj_mat = adj_mat / adj_mat.sum(axis=0)
    id_mat = np.identity(len(nodelist))
    F_mat = 0.85*sla.inv((id_mat-(1-0.85)*adj_mat))
    print('Time it took to make the influence matrix',time.time()-tic)
    np.save('../data_backend2/Influence/%s_data.npy'%net_type,F_mat)

    tic = time.time()
    data = np.load('../data_backend2/Adjacency/%s_data.npy'%net_type)
    print('Time it took to load adj_mat',time.time()-tic)
    
    total_elements = data.shape[0]*data.shape[1]
    num_zeros = total_elements - np.count_nonzero(data)
    frac_zeros = num_zeros/total_elements
    print('For AdjMat',total_elements,num_zeros,frac_zeros)

    tic = time.time()
    data = np.load('../data_backend2/Influence/%s_data.npy'%net_type)
    print('Time it took to load influence matrix',time.time()-tic)

    total_elements = data.shape[0]*data.shape[1]
    num_zeros = total_elements - np.count_nonzero(data)
    frac_zeros = num_zeros/total_elements
    print('For Influence',total_elements,num_zeros,frac_zeros)
    print()