"""Loaders for various data used by GenePlexus."""
import os.path as osp
from typing import Literal

import numpy as np

from . import config
from . import util


def _load_file(
    file_loc: str,
    file_name: str,
    load_method: Literal["npy", "txt"],
) -> np.ndarray:
    file_path = osp.join(file_loc, file_name)
    util.check_file(file_path)

    if load_method == "npy":
        return np.load(file_path)
    elif load_method == "txt":
        return np.loadtxt(file_path, dtype=str)
    else:
        raise ValueError(f"Unknwon load method: {load_method!r}")


def load_node_order(file_loc: str, net_type: config.NET_TYPE) -> np.ndarray:
    """Load network genes."""
    file_name = f"NodeOrder_{net_type}.txt"
    return _load_file(file_loc, file_name, load_method="txt")


def load_genes_universe(
    file_loc: str,
    GSC: config.GSC_TYPE,
    net_type: config.NET_TYPE,
) -> np.ndarray:
    """Load gene universe a given network and GSC."""
    file_name = f"GSC_{GSC}_{net_type}_universe.txt"
    return _load_file(file_loc, file_name, load_method="txt")


def load_gene_features(
    file_loc: str,
    features: config.FEATURE_TYPE,
    net_type: config.NET_TYPE,
) -> np.ndarray:
    """Load gene features."""
    file_name = f"Data_{features}_{net_type}.npy"
    return _load_file(file_loc, file_name, load_method="npy")


def load_correction_order(
    file_loc: str,
    target_set: config.GSC_TYPE,
    net_type: config.NET_TYPE,
) -> np.ndarray:
    """Load correction matrix order."""
    file_name = f"CorrectionMatrixOrder_{target_set}_{net_type}.txt"
    return _load_file(file_loc, file_name, load_method="txt")


def load_correction_mat(
    file_loc: str,
    GSC: config.GSC_TYPE,
    target_set: config.GSC_TYPE,
    net_type: config.NET_TYPE,
    features: config.FEATURE_TYPE,
) -> np.ndarray:
    """Load correction matrix."""
    file_name = f"CorrectionMatrix_{GSC}_{target_set}_{net_type}_{features}.npy"
    return _load_file(file_loc, file_name, load_method="npy")
