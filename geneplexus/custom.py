"""Helper functions for setting up custom networks and GSCs."""
import json
import os.path as osp

import numpy as np

from ._config import logger


def edgelist_to_nodeorder(
    edgelist_loc: str,
    data_dir: str,
    net_name: str,
    sep: str = "\t",
    skiplines: int = 0,
):
    """Convert :term:`edgelist` to node order.

    Args:
        edgelist_loc: Location of the edgelist
        data_dir: The directory to save the file
        net_name: The name of the network
        sep: The separation used in the edgelist file (default tab)
        skiplines: The number of lines to skip for header

    """
    logger.info("Making the NodeOrder File")
    with open(edgelist_loc, "r") as f:
        nodeset = set()
        for idx, line in enumerate(f):
            if idx - skiplines < 0:
                continue
            else:
                nodeset.update(line.strip().split(sep)[:2])
    outfile = osp.join(data_dir, f"NodeOrder_{net_name}.txt")
    logger.info(f"Saving NodeOrder file to {outfile}")
    np.savetxt(outfile, sorted(nodeset), fmt="%s")


def edgelist_to_matrix(
    edgelist_loc: str,
    data_dir: str,
    net_name: str,
    features: str,
    beta: float = 0.85,
    sep: str = "\t",
    skiplines: int = 0,
):
    """Convert :term:`edgelist` to adjacency matrix.

    Args:
        edgelist_loc: Location of the edgelist
        data_dir: The directory to save the file
        net_name: The name of the network
        features: Features for the networks (Adjacency or Influence, All)
        beta: Restart parameter.
        sep: The separation used in the edgelist file (default tab)
        skiplines: The number of lines to skip for header

    """
    if beta < 0 or beta > 1:
        raise ValueError(f"Restart parameter (beta) must be between 0 and 1, got {beta!r}")

    # Load in the NodeOrder file and make node index map
    nodeorder_loc = osp.join(data_dir, f"NodeOrder_{net_name}.txt")
    nodelist = np.loadtxt(nodeorder_loc, dtype=str)
    node_to_ind = {j: i for i, j in enumerate(nodelist)}

    # Make adjacency matrix
    logger.info("Making the adjacency matrix")
    adj_mat = np.zeros((len(nodelist), len(nodelist)), dtype=float)
    with open(edgelist_loc, "r") as f:
        for idx, line in enumerate(f):
            if idx - skiplines < 0:
                continue
            terms = line.strip().split(sep)
            node1, node2 = terms[:2]
            if len(terms) > 3:
                raise ValueError("Too many columns in edgelist file")
            if (node1 not in node_to_ind) or (node2 not in node_to_ind):
                raise KeyError(f"Nodes in Edgelist but not in NodeOrder file ({node1!r} or {node2!r})")
            i, j = node_to_ind[node1], node_to_ind[node2]
            weight = 1.0 if len(terms) == 2 else terms[2]
            adj_mat[i, j] = adj_mat[j, i] = weight

    # Optionally make influence matrix
    if (features == "Influence") or (features == "All"):
        logger.info("Making the influence matrix")
        adj_mat_norm = adj_mat / adj_mat.sum(axis=0)
        id_mat = np.identity(len(nodelist))
        F_mat = beta * np.linalg.inv(id_mat - (1 - beta) * adj_mat_norm)

    # Save the data
    logger.info("Saving the data")
    if (features == "Adjacency") or (features == "All"):
        np.save(osp.join(data_dir, f"Data_Adjacency_{net_name}.npy"), adj_mat)
    if (features == "Influence") or (features == "All"):
        np.save(osp.join(data_dir, f"Data_Influence_{net_name}.npy"), F_mat)


def subset_gsc_to_network(
    data_dir: str,
    net_name: str,
    gsc_name: str,
    max_size: int = 200,
    min_size: int = 10,
):
    """Subset :term:`GSC` using network genes.

    Note:
        Use the :meth:`geneplexus.download.download_select_data` function to
        get the preprocessed GO and DisGeNet files first.

        The NodeOrder file needs to be a single column text file. If not
        supplying custom GSC, the file needs to be in Entrez ID space.

    Args:
        data_dir: The directory to save the file
        net_name: The name of the network
        gsc_name: The name of the GSC
        max_size: Maximum geneset size.
        max_size: Minimum geneset size.

    """
    logger.info("Subsetting the GSC (this make take a few minutes)")
    # load in the NodeOrder file
    nodeorder_loc = osp.join(data_dir, f"NodeOrder_{net_name}.txt")
    nodelist = np.loadtxt(nodeorder_loc, dtype=str)
    # load the orginal GSC
    with open(osp.join(data_dir, f"GSCOriginal_{gsc_name}.json"), "r") as handle:
        gsc_orig = json.load(handle)
    # subset GSc based on network
    universe_genes = np.array([])
    gsc_subset = {}
    for akey in gsc_orig:
        org_genes = gsc_orig[akey]["Genes"]
        genes_tmp = np.intersect1d(nodelist, org_genes)
        if (len(genes_tmp) <= max_size) and (len(genes_tmp) >= min_size):
            gsc_subset[akey] = {"Name": gsc_orig[akey]["Name"], "Genes": genes_tmp.tolist()}
            universe_genes = np.union1d(universe_genes, genes_tmp)
    logger.info("Saving the data")
    with open(osp.join(data_dir, f"GSC_{gsc_name}_{net_name}_GoodSets.json"), "w") as f:
        json.dump(gsc_subset, f, ensure_ascii=False, indent=4)
    np.savetxt(osp.join(data_dir, f"GSC_{gsc_name}_{net_name}_universe.txt"), universe_genes, fmt="%s")
