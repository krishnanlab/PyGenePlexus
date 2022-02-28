import os.path as osp
import pickle
from typing import Dict
from typing import List
from typing import Literal
from typing import Union

import numpy as np

from . import config


def check_file(path: str):
    """Check existence of a file.

    Args:
        path (str): Path to the file.

    Raise:
        FileNotFoundError: if file not exist.

    """
    if not osp.isfile(path):
        raise FileNotFoundError(path)


def read_gene_list(
    path: str,
    sep: str = ", ",
) -> List[str]:
    """Read gene list from flie.

    Args:
        path (str): Path to the input gene list file.
        sep (str): Seperator between genes.

    """
    return [gene.strip("'") for gene in open(path, "r").read().split(sep)]


def load_geneid_conversion(
    file_loc: str,
    src_id_type: config.ID_SRC_TYPE,
    dst_id_type: config.ID_DST_TYPE,
    upper: bool = False,
) -> Dict[str, List[str]]:
    """Load the gene ID conversion mapping.

    Args:
        file_loc (str): Directory containig the ID conversion file.
        src_id_type (ID_SRC_TYPE): Souce gene ID type.
        dst_id_type (ID_DST_TYPE): Destination gene ID type.
        upper (bool): If set to True, then convert all keys to upper case.

    """
    if (src_id_type, dst_id_type) not in config.VALID_ID_CONVERSION:
        raise ValueError(f"Invalid ID conversion from {src_id_type} to {dst_id_type}")

    file_name = f"IDconversion_Homo-sapiens_{src_id_type}-to-{dst_id_type}.pickle"
    file_path = osp.join(file_loc, file_name)
    check_file(file_path)

    with open(file_path, "rb") as handle:
        conversion_map = pickle.load(handle)

    if upper:
        conversion_map = {src.upper(): dst for src, dst in conversion_map.items()}

    return conversion_map


def load_gsc(
    file_loc: str,
    GSC: config.GSC_TYPE,
    net_type: config.NET_TYPE,
) -> Dict[str, Dict[Literal["Name", "Genes"], Union[str, np.ndarray]]]:
    file_name = f"GSC_{GSC}_{net_type}_GoodSets.pickle"
    file_path = osp.join(file_loc, file_name)
    check_file(file_path)

    with open(file_path, "rb") as handle:
        good_sets = pickle.load(handle)

    return good_sets


def load_pretrained_weights(
    file_loc: str,
    target_set: config.GSC_TYPE,
    net_type: config.NET_TYPE,
    features: config.FEATURE_TYPE,
) -> Dict[str, Dict[Literal["Name", "Weights", "PosGenes"], Union[str, np.ndarray]]]:
    file_name = f"PreTrainedWeights_{target_set}_{net_type}_{features}.pickle"
    file_path = osp.join(file_loc, file_name)
    check_file(file_path)

    with open(file_path, "rb") as handle:
        weights_dict = pickle.load(handle)

    return weights_dict


def _load_np_file(
    file_loc: str,
    file_name: str,
    load_method: Literal["npy", "txt"],
) -> np.ndarray:
    """Check np file existence and load."""
    file_path = osp.join(file_loc, file_name)
    check_file(file_path)

    if load_method == "npy":
        return np.load(file_path)
    elif load_method == "txt":
        return np.loadtxt(file_path, dtype=str)
    else:
        raise ValueError(f"Unknwon load method: {load_method!r}")


def load_node_order(file_loc: str, net_type: config.NET_TYPE) -> np.ndarray:
    """Load network genes."""
    file_name = f"NodeOrder_{net_type}.txt"
    return _load_np_file(file_loc, file_name, load_method="txt")


def load_genes_universe(
    file_loc: str,
    GSC: config.GSC_TYPE,
    net_type: config.NET_TYPE,
) -> np.ndarray:
    """Load gene universe a given network and GSC."""
    file_name = f"GSC_{GSC}_{net_type}_universe.txt"
    return _load_np_file(file_loc, file_name, load_method="txt")


def load_gene_features(
    file_loc: str,
    features: config.FEATURE_TYPE,
    net_type: config.NET_TYPE,
) -> np.ndarray:
    """Load gene features."""
    file_name = f"Data_{features}_{net_type}.npy"
    return _load_np_file(file_loc, file_name, load_method="npy")


def load_correction_order(
    file_loc: str,
    target_set: config.GSC_TYPE,
    net_type: config.NET_TYPE,
) -> np.ndarray:
    """Load correction matrix order."""
    file_name = f"CorrectionMatrixOrder_{target_set}_{net_type}.txt"
    return _load_np_file(file_loc, file_name, load_method="txt")


def load_correction_mat(
    file_loc: str,
    GSC: config.GSC_TYPE,
    target_set: config.GSC_TYPE,
    net_type: config.NET_TYPE,
    features: config.FEATURE_TYPE,
) -> np.ndarray:
    """Load correction matrix."""
    file_name = f"CorrectionMatrix_{GSC}_{target_set}_{net_type}_{features}.npy"
    return _load_np_file(file_loc, file_name, load_method="npy")
