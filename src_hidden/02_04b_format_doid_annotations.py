import numpy as np
import pandas as pd
import pickle

fp_hid = '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/data_hidden2/GSCs/'

doid_name_dict = {}            
with open(fp_hid + 'doid_terms.tsv','r') as f:
    for idx, line in enumerate(f):
        line = line.strip().split('\t')
        doid_name_dict[line[0]] = line[1]
print('The number of terms in terms is', len(doid_name_dict))

CO_DO_dict = {}
with open(fp_hid + 'disease_mappings.tsv','r') as f:
    for idx, line in enumerate(f):
        if idx == 0:
            continue
        line = line.strip().split('\t')
        COID = line[0]
        DisID = line[2]+'ID:'+line[3]
        if DisID in doid_name_dict:
            CO_DO_dict[COID] = DisID
print('The len of CO_DO_dict is',len(CO_DO_dict))

results = []
with open(fp_hid + 'curated_gene_disease_associations.tsv','r') as f:
    for idx, line in enumerate(f):
        if idx == 0:
            continue
        line = line.strip().split('\t')
        geneID = line[0]
        diseaseID = line[4]
        if diseaseID not in CO_DO_dict:
            continue
        DOID = CO_DO_dict[diseaseID]
        if DOID not in doid_name_dict:
            continue
        results.append([geneID,DOID,doid_name_dict[DOID]])
df = pd.DataFrame(results,columns=['EntrezID','DOID','Term'])
print('The shape of df is',df.shape)
df.to_csv(fp_hid+'direct_annotations_doid.tsv',sep='\t',header=True,index=False)


