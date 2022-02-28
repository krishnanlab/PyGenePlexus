"""Global variables used by the GenePlexus library."""
from typing import Literal
from typing import Set
from typing import Tuple

ALL_TASKS = ["IDconversion", "MachineLearning", "Similarities", "NetworkGraph"]
ALL_NETWORKS = ["BioGRID", "STRING", "STRING-EXP", "GIANT-TN"]
ALL_FEATURES = ["Adjacency", "Embedding", "Influence"]
ALL_GSCS = ["GO", "DisGeNet"]

URL_AZURE = "https://mancusogeneplexusstorage.blob.core.windows.net/mancusoblob2/"

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
