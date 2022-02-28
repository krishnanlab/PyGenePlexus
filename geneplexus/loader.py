"""Loaders for various data used by GenePlexus."""
import os.path as osp

import numpy as np

from . import config


def load_node_order(file_loc: str, net_type: config.NET_TYPE) -> np.ndarray:
    """Load network genes."""
    file_path = osp.join(file_loc, f"NodeOrder_{net_type}.txt")
    node_order = np.loadtxt(file_path, dtype=str)
    return node_order


def load_genes_universe(
    file_loc: str,
    GSC: config.GSC_TYPE,
    net_type: config.NET_TYPE,
) -> np.ndarray:
    """Load gene universe a given network and GSC."""
    file_path = osp.join(file_loc, f"GSC_{GSC}_{net_type}_universe.txt")
    genes = np.loadtxt(file_path, dtype=str)
    return genes


def load_gene_features(
    file_loc: str,
    features: config.FEATURE_TYPE,
    net_type: config.NET_TYPE,
) -> np.ndarray:
    """Load gene features."""
    file_path = osp.join(file_loc, f"Data_{features}_{net_type}.npy")
    gene_features = np.load(file_path)
    return gene_features


def load_correction_order(
    file_loc: str,
    target_set: config.GSC_TYPE,
    net_type: config.NET_TYPE,
) -> np.ndarray:
    """Load correction matrix order."""
    file_path = osp.join(file_loc, f"CorrectionMatrixOrder_{target_set}_{net_type}.txt")
    correction_order = np.loadtxt(file_path, dtype=str)
    return correction_order


def load_correction_mat(
    file_loc: str,
    GSC: config.GSC_TYPE,
    target_set: config.GSC_TYPE,
    net_type: config.NET_TYPE,
    features: config.FEATURE_TYPE,
) -> np.ndarray:
    """Load correction matrix."""
    file_path = osp.join(file_loc, f"CorrectionMatrix_{GSC}_{target_set}_{net_type}_{features}.npy")
    correction_mat = np.load(file_path)
    return correction_mat
