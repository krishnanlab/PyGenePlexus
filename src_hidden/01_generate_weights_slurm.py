import numpy as np
import pandas as pd
import os
import pickle
import time


tic = time.time()


# dir to put the slurm files
slurm_dir = '/mnt/gs18/scratch/users/mancus16/GenePlexus/slurms/'

combos = []
for anet in ['BioGRID','STRING-EXP','STRING','GIANT-TN']:
    for afeat in ['Embedding','Adjacency','Influence']:
        for agsc in ['DisGeNet','GO']:
            with open('../data_backend/GSCs/%s_%s_GoodSets.pickle'%(agsc,anet),'rb') as handle:
                good_sets = pickle.load(handle)
                if afeat == 'Embedding':
                    set_size = 575
                else:
                    set_size = 48
                num_bins = int(np.ceil(len(good_sets)/set_size))
                for idx in range(num_bins):
                    start_ind = idx*set_size
                    stop_ind = idx*set_size + set_size
                    if stop_ind > len(good_sets):
                        stop_ind = len(good_sets)
                    combos.append([anet,afeat,agsc,start_ind,stop_ind])
print('the number of jobs to submit is',len(combos))              

for idx, acombo in enumerate(combos):

    mylist = ['#!/bin/bash']
    mylist.append('### define resources needed:')
    mylist.append('#SBATCH --time=03:50:00')
    mylist.append('#SBATCH --nodes=1')
    mylist.append('#SBATCH --mem=50G')
    mylist.append('#SBATCH --cpus-per-task=1')
    mylist.append('#SBATCH --job-name=%i-Weight'%idx)
    mylist.append('#SBATCH --account=wang-krishnan')
    mylist.append('#SBATCH --output=%sslurm-%%x-%%j.out'%slurm_dir)
    mylist.append('#SBATCH --exclude=qml-[000-005]')
    mylist.append('module load powertools')
    mylist.append('umask g+rw')
    mylist.append('module use /mnt/research/compbio/krishnanlab/modules/')
    mylist.append('export PATH="/mnt/home/mancus16/software/anaconda3/bin:$PATH"')
    mylist.append('cd /mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/src_hidden')
    mylist.append('python generate_weights.py -n %s -f %s -g %s -stm %i -spm %i'%(acombo[0],acombo[1],acombo[2],acombo[3],acombo[4]))

    with open(slurm_dir + 'weights-%i.sb'%idx, 'w') as thefile:
        for item in mylist:
            thefile.write("%s\n" % item)

    os.system('sbatch ' + slurm_dir + 'weights-%i.sb'%idx)


print('This script took %i minutes to run '%((time.time()-tic)/60))
