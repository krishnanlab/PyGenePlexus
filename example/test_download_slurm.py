import datetime
import glob
import os
import time
from subprocess import PIPE
from subprocess import Popen

import numpy as np

tic = time.time()

myloc = "Azure"
myamount = "full"
num_runs = 5

jobs_final = []
for i in range(num_runs):
    jobs_final.append([myloc, myamount, f"dir{i}"])

# dir to put the slurm files
slurm_dir = f"/mnt/home/mancus16/GenePlexus/speedtest/zipped/{myloc}/"

for idx, ajob in enumerate(jobs_final):
    # make the script
    mylist = ["#!/bin/bash"]
    mylist.append("### define resources needed:")
    mylist.append("#SBATCH --time=03:50:00")
    mylist.append("#SBATCH --nodes=1")
    mylist.append("#SBATCH --cpus-per-task=1")
    mylist.append("#SBATCH --mem=50G")
    mylist.append("#SBATCH --job-name=%s-%s-%s" % (ajob[0], ajob[1], ajob[2]))
    mylist.append("#SBATCH --account=mancuso")
    mylist.append("#SBATCH --output=%sslurm-%%x-%%j.out" % slurm_dir)
    mylist.append("cd")
    mylist.append("module load powertools")
    mylist.append("umask g+rw")
    mylist.append("module use /mnt/research/compbio/krishnanlab/modules/")
    mylist.append("source .bashrc")
    mylist.append('export PATH="/mnt/home/mancus16/software/anaconda3/bin:$PATH"')
    mylist.append("conda activate pygpdpwnloadtest")
    mylist.append("cd /mnt/home/mancus16/GenePlexus/PyGeneplexus/example")
    mylist.append("which python")
    mylist.append(f"python test_download.py -data_loc {ajob[0]} -amount {ajob[1]} -dirname {ajob[2]}")

    with open(slurm_dir + "%s-%s-%s.sb" % (ajob[0], ajob[1], ajob[2]), "w") as thefile:
        for item in mylist:
            thefile.write("%s\n" % item)

    os.system("sbatch " + slurm_dir + "%s-%s-%s.sb" % (ajob[0], ajob[1], ajob[2]))

    p1 = Popen(["squeue", "-u", "mancus16"], stdout=PIPE)
    p2 = Popen(["wc", "-l"], stdin=p1.stdout, stdout=PIPE)
    # njobs has like 7 extra lines becuase of the header, but it is good to over estimate the number of jobs in the queue
    njobs, err = p2.communicate()
    njobs = int(njobs)

    while njobs > 800:
        print("More than 800 jobs in queue")
        time.sleep(360)
        p1 = Popen(["squeue", "-u", "mancus16"], stdout=PIPE)
        p2 = Popen(["wc", "-l"], stdin=p1.stdout, stdout=PIPE)
        njobs, err = p2.communicate()
        njobs = int(njobs)

print("This script took %i minutes to run " % ((time.time() - tic) / 60))
