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
SPECIES_TYPE = Literal["Human", "Mouse", "Fly", "Worm", "Fish", "Yeast"]

NET_SELECTION_TYPE = Union[Literal["All"], NET_TYPE, List[NET_TYPE]]
FEATURE_SELECTION_TYPE = Union[Literal["All"], FEATURE_TYPE, List[FEATURE_TYPE]]
GSC_SELECTION_TYPE = Union[Literal["All"], GSC_TYPE, List[GSC_TYPE]]
SPECIES_SELECTION_TYPE = Union[Literal["All"], SPECIES_TYPE, List[SPECIES_TYPE]]

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
]
