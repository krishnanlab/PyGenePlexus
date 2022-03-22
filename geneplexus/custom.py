import json
import os
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
    """Convert edge list to node order.

    Note:
        The edgelist file needs to be two or three columns. The first two
        columns being the edges. The third column is assumed to be the edge
        weight if it exsits. If not supplying custom GSC, the file needs to be
        in Entrez ID space.

    Args:
        edgelist_loc: Location of the edgelist
        data_dir: The directory to save the file
        net_name: The name of the network
        sep: The separation used in the edgelist file (default tab)
        skiplines: The number of lines to skip for header

    """
    logger.info("Making the NodeOrder File")
    with open(edgelist_loc, "r") as f:
        nodelist = set()
        for idx, line in enumerate(f):
            if idx - skiplines < 0:
                continue
            else:
                line = line.strip().split(sep)
                nodelist.add(line[0])
                nodelist.add(line[1])
    outfile = osp.join(data_dir, f"NodeOrder_{net_name}.txt")
    logger.info(f"Saving NodeOrder file to {outfile}")
    np.savetxt(outfile, list(nodelist), fmt="%s")


def edgelist_to_matrix(
    edgelist_loc: str,
    nodeorder_loc: str,
    data_dir: str,
    net_name: str,
    features: str,
    sep: str = "\t",
    skiplines: int = 0,
):
    """Convert edge list to adjacency matrix.

    Note:
        The edgelist file needs to be two or three columns. The first two
        columns being the edges. The third column is assumed to be the edge
        weight if it exsits. If not supplying custom GSC, the file needs to be
        in Entrez ID space. Finally, the NodeOrder file needs to be a single
        column text file.

    Args:
        edgelist_loc: Location of the edgelist
        nodeorder_loc: Location of the NodeOrder file
        data_dir: The directory to save the file
        net_name: The name of the network
        features: Features for the networks (Adjacency or Influence, All)
        sep: The separation used in the edgelist file (default tab)
        skiplines: The number of lines to skip for header

    """
    # load in the NodeOrder file
    nodelist = np.loadtxt(nodeorder_loc, dtype=str)
    # make node to ind dict
    node_to_ind = {}
    for idx, anode in enumerate(nodelist):
        node_to_ind[anode] = idx
    # make adjacency matrix
    logger.info("Making the adjacency matrix")
    adj_mat = np.zeros((len(nodelist), len(nodelist)), dtype=float)
    with open(edgelist_loc, "r") as f:
        for idx, line in enumerate(f):
            if idx - skiplines < 0:
                continue
            line = line.strip().split(sep)
            if (line[0] not in node_to_ind) or (line[1] not in node_to_ind):
                raise KeyError("Nodes in Edgelist not in NodeOrder file")
            if len(line) == 2:
                adj_mat[node_to_ind[line[0]], node_to_ind[line[1]]] = 1.0
                adj_mat[node_to_ind[line[1]], node_to_ind[line[0]]] = 1.0
            elif len(line) == 3:
                adj_mat[node_to_ind[line[0]], node_to_ind[line[1]]] = float(line[2])
                adj_mat[node_to_ind[line[1]], node_to_ind[line[0]]] = float(line[2])
            else:
                raise ValueError("Too many columns in edgelist file")
    if (features == "Influence") or (features == "All"):
        logger.info("Making the influence matrix")
        # make influence matrix
        adj_mat_norm = adj_mat / adj_mat.sum(axis=0)
        id_mat = np.identity(len(nodelist))
        F_mat = 0.85 * sla.inv(id_mat - (1 - 0.85) * adj_mat_norm)
    # save the data
    logger.info("Saving the data")
    if (features == "Adjacency") or (features == "All"):
        np.save(osp.join(data_dir, f"Data_Adjacency_{net_name}.npy"), adj_mat)
    if (features == "Influence") or (features == "All"):
        np.save(osp.join(data_dir, f"Data_Influence_{net_name}.npy"), F_mat)


def subset_GSC_to_network(
    nodeorder_loc: str,
    data_dir: str,
    GSC_name: str,
):
    """Subset geneset collection using network genes.

    Note:
        Use the :meth:`geneplexus.download.download_select_data` function to
        get the preprocessed GO and DisGeNet files first.

        The NodeOrder file needs to be a single column text file. If not
        supplying custom GSC, the file needs to be in Entrez ID space.

    Args:
        nodeorder_loc: Location of the NodeOrder file
        data_dir: The directory to save the file
        GSC_name: The name of the GSC

    """
    logger.info("Subsetting the GSC")
    logger.info("This make take a few minutes")
    # load in the NodeOrder file
    nodelist = np.loadtxt(nodeorder_loc, dtype=str)
    # load the orginal GSC
    with open(osp.join(data_dir, f"GSCOriginal_{GSC_name}.json"), "r") as handle:
        GSCorg = json.load(handle)
    # subset GSc based on network
    universe_genes = np.array([])
    GSCsubset = {}
    for akey in GSCorg:
        org_genes = GSCorg[akey]["Genes"]
        genes_tmp = np.intersect1d(nodelist, org_genes)
        if (len(genes_tmp) <= 200) and (len(genes_tmp) >= 10):
            GSCsubset[akey] = {"Name": GSCorg[akey]["Name"], "Genes": genes_tmp.tolist()}
            universe_genes = np.union1d(universe_genes, genes_tmp)
    logger.info("Saving the data")
    net_name = os.path.basename(nodeorder_loc).split("_")[1].split(".t")[0]
    with open(osp.join(data_dir, f"GSC_{GSC_name}_{net_name}_GoodSets.json"), "w") as f:
        json.dump(GSCsubset, f, ensure_ascii=False, indent=4)
    np.savetxt(osp.join(data_dir, f"GSC_{GSC_name}_{net_name}_universe.txt"), universe_genes, fmt="%s")
