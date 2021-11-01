import numpy as np
import pickle
import os
import time
import glob
import datetime
import sys

'''
I noticed at least one case where a gene in the zebrafish BioGRID data was actually a
human gene. I'm going to only inlcude genes that we have MyGeneInfo for in that species.
This takes extra time but is worth it

I also tried doing both ALL and ORGANISM-Home-sapaiens and the end result is the same

'''
    
def check_exists(FileNameList):
    for item in FileNameList:
        if os.path.isfile(item) == True:
            raise FileExistsError('{} already exists'.format(item))
            
def get_MyGeneInfo_IDs():
    '''
    This function will get the Entrez IDs we got from NCBI
    '''
    fp_MGI = '/mnt/research/compbio/krishnanlab/data/MyGeneInfo/20201029_Entrez_Multiple-Species/Homo_sapiens/'
    FileNames = glob.glob(fp_MGI+'*.json')
    IDs = [int(item.split('/')[-1].split('.')[0].split('-')[1]) for item in FileNames]
    return IDs
            
def main(File_org,IDs):
    Mod_list = []
    Unq_genes = []
    Unq_genes_org = []
    cnt1 = 0
    with open(File_org, mode='r') as f_old:
        for idx, line in enumerate(f_old):
            if idx == 0:
                continue
            else:
                cnt1 = cnt1 + 1
            try:
                EntA = int(line.split('\t')[0].split(':')[1])
                EntB = int(line.split('\t')[1].split(':')[1])
                Unq_genes_org.append(EntA)
                Unq_genes_org.append(EntB)
            except ValueError:
                continue
            if (EntA not in IDs) or (EntB not in IDs):
                continue
            else:
                pass
            if EntB < EntA:
                EntB,EntA = EntA,EntB
            else:
                pass
            if EntA == EntB:
                continue
            else:
                Mod_list.append((EntA,EntB))
                Unq_genes.append(EntA)
                Unq_genes.append(EntB)
    Mod_list = list(set(Mod_list))   
    Unq_genes = list(set(Unq_genes))
    Unq_genes_org = list(set(Unq_genes_org))          
    return Mod_list, Unq_genes, cnt1, Unq_genes_org
                        
def make_edgelist(Mod_list,File_Mod): 
    mylist = sorted(Mod_list, key=lambda element: (element[0], element[1]))
    myarr = np.array(mylist)
    np.savetxt(File_Mod,myarr,fmt=['%i','%i'],delimiter='\t')


if __name__ == "__main__":
    
    # This will make an unweighted edge list as the scores in BioGrid a more complicated
    # to convert into meaniful edgeweights
    # set file paths to read and save networks
    fp = '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/data_hidden2/Network_Downloads/BioGRID/'
    File_org = fp + 'BIOGRID-ORGANISM-Homo_sapiens-4.2.191.mitab.txt'
    # File_org = fp + 'BIOGRID-ALL-4.2.191.mitab.txt'
    fp2 = '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/data_backend2/Edgelists/'
    File_Mod = fp2 + 'BioGRID.edg'
    now = datetime.datetime.now()
    sys.stdout = open(fp+'%04d%02d%02d_runlog.txt'%(now.year,now.month,now.day), 'w')
    # check if the modified files exist
    FileList = [File_Mod]
    # check_exists(FileList)
    # get MyGeneInfo Genes
    IDs = get_MyGeneInfo_IDs()
    print('There are %i genes for Homo_sapiens in MyGeneInfo'%len(IDs))
    # # # run main file
    print('Finding valid edges')
    tic = time.time()
    Mod_list, Unq_genes, cnt1, Unq_genes_org = main(File_org,IDs)
    print('The number of edges in the original network is',cnt1)
    print('The number of unique genes in the original network is', len(Unq_genes_org))
    print('The number of edges in the modified network is', len(Mod_list))
    print('The number of unique genes in the modified network is', len(Unq_genes))
    print('The number of seconds it took to find valid edges is',int(time.time()-tic))
    # make edgelist
    print('Writing final edgelist for')
    tic = time.time()
    make_edgelist(Mod_list,File_Mod)
    print('The number of seconds it took to write the edgelist is',int(time.time()-tic))
    