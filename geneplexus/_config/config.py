"""Global variables used by the GenePlexus library."""
import os.path as osp
import pathlib
from typing import Any
from typing import Dict
from typing import List
from typing import Literal
from typing import Set
from typing import Tuple
from typing import Union

import numpy as np

MAX_RETRY = 10  # maximum number of retries for downloading

URL_DICT = {
    "Zenodo": "https://zenodo.org/records/14750555/files/",
    "ZenodoAPI": "https://zenodo.org/api/records/14750555/files/",
}

CONFIG_PATH = pathlib.Path(__file__).parent.absolute()
DATA_FILENAMES_PATH = osp.join(CONFIG_PATH, "data_filenames.txt")

LOG_LEVELS = ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"]

ALL_NETWORKS = ["BioGRID", "STRING", "IMP"]
ALL_FEATURES = ["SixSpeciesN2V"]
ALL_GSCS = ["GO", "Monarch", "Mondo", "Combined"]
ALL_SPECIES = ["Human", "Mouse", "Fly", "Worm", "Zebrafish", "Yeast"]
ALL_CLUSTERING = ["louvain", "domino"]

DEFAULT_LOGREG_KWARGS: Dict[str, Any] = {
    "max_iter": 10000,
    "solver": "lbfgs",
    "penalty": "l2",
    "C": 1.0,
}

LOG_LEVEL_TYPE = Literal["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"]

ID_SRC_TYPE = Literal["ENSG", "ENSP", "ENST", "Entrez", "Symbol"]
ID_DST_TYPE = Literal["Entrez", "ENSG", "Name", "Symbol"]
VALID_ID_CONVERSION: Set[Tuple[ID_SRC_TYPE, ID_DST_TYPE]] = {
    ("ENSG", "Entrez"),
    ("ENSP", "Entrez"),
    ("ENST", "Entrez"),
    ("Symbol", "Entrez"),
    ("Entrez", "ENSG"),
    ("Entrez", "Name"),
    ("Entrez", "Symbol"),
}

NET_TYPE = Literal["BioGRID", "STRING", "IMP"]
FEATURE_TYPE = Literal["SixSpeciesN2V"]
GSC_TYPE = Literal["GO", "Monarch", "Mondo", "Combined"]
SPECIES_TYPE = Literal["Human", "Mouse", "Fly", "Worm", "Zebrafish", "Yeast"]
CLUSTERING_TYPE = Literal["louvain", "domino"]

NET_SELECTION_TYPE = Union[Literal["All"], NET_TYPE, List[NET_TYPE]]
FEATURE_SELECTION_TYPE = Union[Literal["All"], FEATURE_TYPE, List[FEATURE_TYPE]]
SPECIES_SELECTION_TYPE = Union[Literal["All"], SPECIES_TYPE, List[SPECIES_TYPE]]
GSC_SELECTION_TYPE = Union[GSC_TYPE, List[GSC_TYPE]]  # slightly different than the others
CLUSTERING_SELECTION_TYPE = Union[CLUSTERING_TYPE, List[CLUSTERING_TYPE]]

ID_CONVERSION_MAP_TYPE = Dict[str, List[str]]
GSC_DATA_TYPE = Dict[str, Dict[Literal["Name", "Genes"], Union[str, np.ndarray]]]
BIOMART_DATA_TYPE = Dict[str, str]
PRETRAINED_DATA_TYPE = Dict[str, Dict[Literal["Name", "Weights", "PosGenes"], Union[str, np.ndarray]]]

COMBINED_CONTEXTS: Dict[str, Any] = {
    "Human": "Combined - Gene Set Contexts Used [GO, Monarch, Mondo]",
    "Mouse": "Combined - Gene Set Contexts Used [GO, Monarch]",
    "Zebrafish": "Combined - Gene Set Contexts Used [GO, Monarch]",
    "Worm": "Combined - Gene Set Contexts Used [GO, Monarch]",
    "Yeast": "Combined - Gene Set Contexts Used [GO, Monarch]",
    "Fly": "Combined - Gene Set Contexts Used [GO]",
}

# Note here, boolean arguments need to be changed in CLI manually if changed here
DEFAULT_PARAMETERS = {
    "file_loc": None,
    "net_type": "STRING",
    "features": "SixSpeciesN2V",
    "sp_trn": "Human",
    "sp_res": "Human",
    "gsc_trn": "Combined",
    "gsc_res": "Combined",
    "input_genes": None,
    "input_negatives": None,
    "auto_download": False,
    "log_level": "INFO",
    "log_to_file": False,
    "clust_method": "louvain",
    "clust_min_size": 5,
    "clust_weighted": True,
    "scale": False,
    "min_num_pos": 5,
    "min_num_pos_cv": 15,
    "num_folds": 3,
    "null_val": None,
    "random_state": 0,
    "cross_validate": True,
    "num_nodes": 50,
    "output_dir": None,
    "save_type": "all",
    "zip_output": False,
    "overwrite": False,
}

__all__ = [
    "URL_DICT",
    "CONFIG_PATH",
    "DATA_FILENAMES_PATH",
    "ALL_NETWORKS",
    "ALL_FEATURES",
    "ALL_GSCS",
    "ALL_SPECIES",
    "LOG_LEVEL_TYPE",
    "ID_SRC_TYPE",
    "ID_DST_TYPE",
    "VALID_ID_CONVERSION",
    "NET_TYPE",
    "FEATURE_TYPE",
    "GSC_TYPE",
    "SPECIES_TYPE",
    "NET_SELECTION_TYPE",
    "FEATURE_SELECTION_TYPE",
    "GSC_SELECTION_TYPE",
    "SPECIES_SELECTION_TYPE",
    "ID_CONVERSION_MAP_TYPE",
    "GSC_DATA_TYPE",
    "PRETRAINED_DATA_TYPE",
    "COMBINED_CONTEXTS",
    "DEFAULT_PARAMETERS",
]
