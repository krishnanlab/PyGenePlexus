import numpy as np
import pandas as pd
import pickle

##################### Disgent section ######################################################

# find all genes in disease collection and all sets of a certain size
print('DisGeNet')
for net_type in ['BioGRID', 'STRING-EXP', 'STRING', 'GIANT-TN']:
    net_genes = np.loadtxt('../data_backend/Node_Orders/%s_nodelist.txt'%net_type,dtype=str)
    all_dis_genes = np.array([])
    universe_dis_genes = np.array([])
    good_dis_sets = {}
    with open('../data_hidden/GSCs/disgenet_disease-genes_propagated.gmt','r') as f:
        for idx, line in enumerate(f):
            line = line.strip().split('\t')
            genes_tmp = np.intersect1d(line[2:],net_genes)
            all_dis_genes = np.union1d(all_dis_genes,genes_tmp)
            if (len(genes_tmp) <= 200) and (len(genes_tmp) >= 10):
                good_dis_sets[line[0]] = {'Name':line[1],'Genes':genes_tmp}
                universe_dis_genes = np.union1d(universe_dis_genes,genes_tmp)
    np.savetxt('../data_backend/GSCs/DisGeNet_%s_universe.txt'%net_type,universe_dis_genes,fmt='%s')
    with open('../data_backend/GSCs/DisGeNet_%s_GoodSets.pickle'%net_type,'wb') as handle:
        pickle.dump(good_dis_sets, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print(net_type)
    print(len(all_dis_genes))
    print(len(universe_dis_genes))
    print(len(good_dis_sets))

print()
print('GO')
for net_type in ['BioGRID', 'STRING-EXP', 'STRING', 'GIANT-TN']:
    net_genes = np.loadtxt('../data_backend/Node_Orders/%s_nodelist.txt'%net_type,dtype=str)
    all_GO_genes = np.array([])
    universe_GO_genes = np.array([])
    good_GO_sets = {}
    with open('../data_hidden/GSCs/GO__propagated-annotations__Homo_sapiens__Entrez__BP__EXP_IDA_IPI_IMP_IGI_TAS_IC.tsv','r') as f:
        for idx, line in enumerate(f):
            if idx == 0:
                continue
            line = line.strip().split('\t')
            genes_tmp = line[3].split(', ')
            genes_tmp = np.intersect1d(genes_tmp,net_genes)
            all_GO_genes = np.union1d(all_GO_genes,genes_tmp)
            if (len(genes_tmp) <= 200) and (len(genes_tmp) >= 10):
                good_GO_sets[line[0]] = {'Name':line[1],'Genes':genes_tmp}
                universe_GO_genes = np.union1d(universe_GO_genes,genes_tmp)
    np.savetxt('../data_backend/GSCs/GO_%s_universe.txt'%net_type,universe_GO_genes,fmt='%s')
    with open('../data_backend/GSCs/GO_%s_GoodSets.pickle'%net_type,'wb') as handle:
        pickle.dump(good_GO_sets, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print(net_type)
    print(len(all_GO_genes))
    print(len(universe_GO_genes))
    print(len(good_GO_sets))

