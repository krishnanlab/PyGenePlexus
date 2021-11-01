import numpy as np
import pandas as pd
import pickle

##################### Disgent section ######################################################

# find all genes in disease collection and all sets of a certain size
print('DisGeNet')
for net_type in ['BioGRID', 'STRING-EXP', 'STRING', 'GIANT-TN']:
    with open('../data_backend/GSCs/DisGeNet_%s_GoodSets.pickle'%net_type,'rb') as handle:
        good_dis_sets = pickle.load(handle)
    print(net_type)
    print(len(good_dis_sets))

print()
print('GO')
for net_type in ['BioGRID', 'STRING-EXP', 'STRING', 'GIANT-TN']:
    with open('../data_backend/GSCs/GO_%s_GoodSets.pickle'%net_type,'rb') as handle:
        good_GO_sets = pickle.load(handle)
    print(net_type)
    print(len(good_GO_sets))

