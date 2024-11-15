import os
import os.path as osp
import pathlib
import shutil
import time

import numpy as np

import geneplexus

# Set up data directory
homedir = pathlib.Path(__file__).absolute().parent
datadir = osp.join(homedir, "data")
os.makedirs(datadir, exist_ok=True)


"""
each file is separated by the species. Can select all by using
species = ["Human", "Mouse", "Fly", "Worm", "Zebrafish", "Yeast"]
or
species = "All"

for a subset just include desired species (example for just Mouse and Human)
species = ["Human", "Mouse"]
"""

geneplexus.download.download_select_data(
    datadir,
    species=["Human", "Mouse", "Fly", "Worm", "Zebrafish", "Yeast"],
)
