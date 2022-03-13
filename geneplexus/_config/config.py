"""Global variables used by the GenePlexus library."""
import os.path as osp
import pathlib
from typing import Dict
from typing import List
from typing import Literal
from typing import Set
from typing import Tuple
from typing import Union

import numpy as np

URL_AZURE = "https://pygeneplexusstacct.blob.core.windows.net/geneplexusblob/"
CONFIG_PATH = pathlib.Path(__file__).parent.absolute()
DATA_FILENAMES_PATH = osp.join(CONFIG_PATH, "data_filenames.txt")

ALL_TASKS = ["IDconversion", "MachineLearning", "Similarities", "NetworkGraph"]
ALL_NETWORKS = ["BioGRID", "STRING", "STRING-EXP", "GIANT-TN"]
ALL_FEATURES = ["Adjacency", "Embedding", "Influence"]
ALL_GSCS = ["GO", "DisGeNet"]

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
NET_TYPE = Literal["BioGRID", "STRING", "STRING-EXP", "GIANT-TN"]
GSC_TYPE = Literal["GO", "DisGeNet"]
FEATURE_TYPE = Literal["Adjacency", "Embedding", "Influence"]

ID_CONVERSION_MAP_TYPE = Dict[str, List[str]]
GSC_DATA_TYPE = Dict[str, Dict[Literal["Name", "Genes"], Union[str, np.ndarray]]]
PRETRAINED_DATA_TYPE = Dict[str, Dict[Literal["Name", "Weights", "PosGenes"], Union[str, np.ndarray]]]

__all__ = [
    "URL_AZURE",
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
    "NET_TYPE",
    "GSC_TYPE",
    "FEATURE_TYPE",
    "ID_CONVERSION_MAP_TYPE",
    "GSC_DATA_TYPE",
    "PRETRAINED_DATA_TYPE",
]
