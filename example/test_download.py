import argparse
import os
import os.path as osp
import pathlib
import shutil
import time

import numpy as np

import geneplexus

parser = argparse.ArgumentParser()
parser.add_argument(
    "-dirname",
    default="dir1",
    type=str,
    help="This dir will be end point",
)
parser.add_argument(
    "-data_loc",
    default="Zenodo",
    type=str,
    help="Where the data is stored (Azure or Zenodo)",
)
parser.add_argument(
    "-amount",
    default="full",
    type=str,
    help="either full or subset",
)
args = parser.parse_args()
dirname = args.dirname
data_loc = args.data_loc
amount = args.amount

dirname_full = f"{data_loc}_{amount}_{dirname}"

tic = time.time()

# Set up directories
homedir = pathlib.Path(__file__).absolute().parent
datadir = osp.join(homedir, "data", dirname_full)
tic1 = time.time()
if os.path.exists(datadir):
    shutil.rmtree(datadir)
os.makedirs(datadir)
print("The time it took to delete/create the directory was", (time.time() - tic1) / 60)

# Get the data
if amount == "full":
    print(f"Start downloading data and saving to: {datadir}")
    geneplexus.download.download_select_data(
        datadir,
        tasks="All",
        networks="All",
        features="All",
        gscs="All",
        data_loc=data_loc,
    )

elif amount == "subset":
    print(f"Start downloading data and saving to: {datadir}")
    geneplexus.download.download_select_data(
        datadir,
        tasks="All",
        networks="BioGRID",
        features="Embedding",
        gscs="GO",
        data_loc=data_loc,
    )

print("Done downlaoding")
print("The time it took in minutes was", (time.time() - tic) / 60)
