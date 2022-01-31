import numpy as np
import pandas as pd
import glob
import shutil

# This script will change how data_back end is structured and labeled
# note that this moves the files. If you want to keep the origianals
# make a copy first

fp_old = '/Users/christophermancuso/Documents/DataSets/data_backend2/'
fp_new = '/Users/christophermancuso/Documents/DataSets/Geneplexus_data/'

# change the adjacency matrix folder
FNs = glob.glob(fp_old+'Adjacency/*.npy')
for aFN in FNs:
    network = aFN.strip().split('/')[-1].split('_')[0]
    new_name = 'Adjacency_%s.npy'%network
    shutil.move(aFN,fp_new+new_name)

# change the embedding matrix folder
FNs = glob.glob(fp_old+'Embedding/*.npy')
for aFN in FNs:
    network = aFN.strip().split('/')[-1].split('_')[0]
    new_name = 'Embedding_%s.npy'%network
    shutil.move(aFN,fp_new+new_name)

# change the influence matrix folder
FNs = glob.glob(fp_old+'Influence/*.npy')
for aFN in FNs:
    network = aFN.strip().split('/')[-1].split('_')[0]
    new_name = 'Influence_%s.npy'%network
    shutil.move(aFN,fp_new+new_name)

# change the node_order folder
FNs = glob.glob(fp_old+'Node_Orders/*.txt')
for aFN in FNs:
    network = aFN.strip().split('/')[-1].split('_')[0]
    new_name = 'NodeOrders_%s.txt'%network
    shutil.move(aFN,fp_new+new_name)

# change the Edgelist folder
FNs = glob.glob(fp_old+'Edgelists/*.edg')
for aFN in FNs:
    network = aFN.strip().split('/')[-1].split('.e')[0]
    new_name = 'Edgelist_%s.edg'%network
    shutil.move(aFN,fp_new+new_name)

# change the Id_conversion folder
FNs = glob.glob(fp_old+'ID_conversion/*.pickle')
for aFN in FNs:
    specie = aFN.strip().split('/')[-1].split('__')[0]
    specie = specie.replace('_','-')
    IDtypes = aFN.strip().split('/')[-1].split('__')[1]
    new_name = 'IDconversion_%s_%s.pickle'%(specie,IDtypes)
    shutil.move(aFN,fp_new+new_name)

# change the GSCs folder
FNs = glob.glob(fp_old+'GSCs/*')
for aFN in FNs:
    old_name = aFN.strip().split('/')[-1]
    new_name = 'GSC_%s'%old_name
    shutil.move(aFN,fp_new+new_name)

# change the PreTrainedModels folder
FNs = glob.glob(fp_old+'PreTrainedModels/*pickle')
for aFN in FNs:
    old_name = aFN.strip().split('/')[-1].split('_Model')[0]
    new_name = 'PreTrainedWeights_%s.pickle'%old_name
    shutil.move(aFN,fp_new+new_name)

# change the CorMats in CorrectionMatrices folder
FNs = glob.glob(fp_old+'CorrectionMatrices/*npy')
for aFN in FNs:
    old_name = aFN.strip().split('/')[-1].split('_CorMat')[0]
    new_name = 'CorrectionMatrix_%s.npy'%old_name
    shutil.move(aFN,fp_new+new_name)

# change the orders in CorrectionMatrices folder
FNs = glob.glob(fp_old+'CorrectionMatrices/*txt')
for aFN in FNs:
    old_name = aFN.strip().split('/')[-1].split('_Orders')[0]
    new_name = 'CorrectionMatrixOrder_%s.txt'%old_name
    shutil.move(aFN,fp_new+new_name)