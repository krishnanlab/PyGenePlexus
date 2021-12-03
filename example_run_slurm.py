import numpy as np
import pandas as pd
import os
import pickle
import time
import glob

# these can be changed
net_type = 'STRING'
features = 'Embedding'
GSC = 'DisGeNet'

# dir to put the slurm files
slurm_dir = 'slurms/'
if not os.path.exists(slurm_dir):
    os.makedirs(slurm_dir)

# get all the files and then get the part we care about for file name
files_to_do = glob.glob('/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/disease_gene_files/*.txt')
jobnames = [item.strip().split('/')[-1].split('.t')[0]+'--%s--%s--%s'%(net_type,features,GSC) for item in files_to_do]


for idx, afile in enumerate(files_to_do):

    mylist = ['#!/bin/bash']
    mylist.append('### define resources needed:')
    mylist.append('#SBATCH --time=03:50:00')
    mylist.append('#SBATCH --nodes=1')
    mylist.append('#SBATCH --mem=50G')
    mylist.append('#SBATCH --cpus-per-task=1')
    mylist.append('#SBATCH --job-name=%s'%jobnames[idx])
    mylist.append('#SBATCH --account=wang-krishnan')
    mylist.append('#SBATCH --output=%sslurm-%%x-%%j.out'%slurm_dir)
    mylist.append('module load powertools')
    mylist.append('umask g+rw')
    mylist.append('export PATH="/mnt/home/mancus16/software/anaconda3/bin:$PATH"') #### change this
    mylist.append('cd /mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GeneplexusPublic') ### change this
    mylist.append('python example_run.py -i %s -j %s -n %s -f %s -g %s'%(afile,jobnames[idx],net_type,features,GSC))

    with open(slurm_dir + 'Geneplexus-%s.sb'%jobnames[idx], 'w') as thefile:
        for item in mylist:
            thefile.write("%s\n" % item)

    os.system('sbatch ' + slurm_dir + 'Geneplexus-%s.sb'%jobnames[idx])
