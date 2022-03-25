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

URL_DATA = "https://zenodo.org/record/6383205/files/"
CONFIG_PATH = pathlib.Path(__file__).parent.absolute()
DATA_FILENAMES_PATH = osp.join(CONFIG_PATH, "data_filenames.txt")

ALL_TASKS = ["IDconversion", "MachineLearning", "Similarities", "NetworkGraph", "OriginalGSCs"]
ALL_NETWORKS = ["BioGRID", "STRING", "STRING-EXP", "GIANT-TN"]
ALL_FEATURES = ["Adjacency", "Embedding", "Influence"]
ALL_GSCS = ["GO", "DisGeNet"]

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
    ("Entrez", "ENSG"),
    ("Entrez", "Name"),
    ("Entrez", "Symbol"),
    ("Symbol", "Entrez"),
}

TASK_TYPE = Literal["IDconversion", "MachineLearning", "Similarities", "NetworkGraph", "OriginalGSCs"]
NET_TYPE = Literal["BioGRID", "STRING", "STRING-EXP", "GIANT-TN"]
FEATURE_TYPE = Literal["Adjacency", "Embedding", "Influence"]
GSC_TYPE = Literal["GO", "DisGeNet"]

TASK_SELECTION_TYPE = Union[Literal["All"], List[TASK_TYPE]]
NET_SELECTION_TYPE = Union[Literal["All"], List[NET_TYPE]]
FEATURE_SELECTION_TYPE = Union[Literal["All"], List[FEATURE_TYPE]]
GSC_SELECTION_TYPE = Union[Literal["All"], List[GSC_TYPE]]

ID_CONVERSION_MAP_TYPE = Dict[str, List[str]]
GSC_DATA_TYPE = Dict[str, Dict[Literal["Name", "Genes"], Union[str, np.ndarray]]]
PRETRAINED_DATA_TYPE = Dict[str, Dict[Literal["Name", "Weights", "PosGenes"], Union[str, np.ndarray]]]

__all__ = [
    "URL_DATA",
    "CONFIG_PATH",
    "DATA_FILENAMES_PATH",
    "ALL_TASKS",
    "ALL_NETWORKS",
    "ALL_FEATURES",
    "ALL_GSCS",
    "LOG_LEVEL_TYPE",
    "ID_SRC_TYPE",
    "ID_DST_TYPE",
    "VALID_ID_CONVERSION",
    "TASK_TYPE",
    "NET_TYPE",
    "FEATURE_TYPE",
    "GSC_TYPE",
    "TASK_SELECTION_TYPE",
    "NET_SELECTION_TYPE",
    "FEATURE_SELECTION_TYPE",
    "GSC_SELECTION_TYPE",
    "ID_CONVERSION_MAP_TYPE",
    "GSC_DATA_TYPE",
    "PRETRAINED_DATA_TYPE",
]
