import copy
import sys

import networkx as nx
import numpy as np
import pcst_fast
from networkx.algorithms.community.centrality import girvan_newman
from networkx.algorithms.community.quality import modularity


def run_pcst(G_cc, n_steps, nodes, prize_factor):
    ## set prize ##
    prizes = get_pcst_prize(G_cc, prize_factor, n_steps)
    vertices_prizes = []
    for cur_node in nodes:
        vertices_prizes.append(
            G_cc.nodes[cur_node]["Active"] if G_cc.nodes[cur_node]["Active"] else prizes[cur_node],
        )
    ## set cost ##
    edges_grid = []
    for cur_edge in G_cc.edges:
        edges_grid.append([nodes.index(cur_edge[0]), nodes.index(cur_edge[1])])
    edges_costs = []
    for cur_edge in edges_grid:
        u_score = 0 if G_cc.nodes[nodes[cur_edge[0]]]["Active"] else 0.9999
        v_score = 0 if G_cc.nodes[nodes[cur_edge[1]]]["Active"] else 0.9999

        edges_costs.append(np.min([u_score, v_score]))
    ## find pcst component by running pcst fast##
    root = -1
    num_clusters = 1
    pruning = "strong"  # 'none'
    verbosity_level = 0
    vertices, edges = pcst_fast.pcst_fast(
        edges_grid,
        vertices_prizes,
        edges_costs,
        root,
        num_clusters,
        pruning,
        verbosity_level,
    )

    return edges, edges_grid


def get_pcst_prize(G_cc, prize_factor, n_steps):
    prizes = {}
    p_cc = linear_threshold(G_cc, [n for n in G_cc.nodes if G_cc.nodes[n]["Active"] > 0], steps=n_steps)
    for p_node in G_cc.nodes:
        prizes[p_node] = 0
    for i_cur_layer, cur_layer in enumerate(p_cc):
        for cur_node in cur_layer:
            prizes[cur_node] += prize_factor**i_cur_layer

    return prizes


def get_putative_modules(
    G_subslice,
    domino_weight,
    improvement_delta,
    modularity_score_objective,
):
    """"""
    G_optimized = G_subslice.copy()

    # clean subslice from cycles and isolated nodes
    G_optimized.remove_edges_from(list(nx.selfloop_edges(G_optimized)))
    G_optimized.remove_nodes_from(list(nx.isolates(G_optimized)))

    # start so if module size < 100 won't be further sliced
    is_enriched_sublice = len(G_optimized.nodes) < 100
    best_modularity = -1
    # loop executes until subslice is enriched
    while not is_enriched_sublice:
        is_enriched_sublice, best_modularity, G_optimized = split_subslice_into_putative_modules(
            G_optimized,
            improvement_delta,
            modularity_score_objective,
            best_modularity,
            domino_weight,
        )
    G_optimized.remove_nodes_from(list(nx.isolates(G_optimized)))
    # get all cc's from final enriched subslice
    cc_optimized = (
        [] if len(G_optimized.nodes) == 0 else [G_optimized.subgraph(c) for c in nx.connected_components(G_optimized)]
    )

    return cc_optimized


def split_subslice_into_putative_modules(
    G_optimized,
    improvement_delta,
    modularity_score_objective,
    best_modularity,
    domino_weight,
):
    cur_components = [G_optimized.subgraph(c) for c in nx.connected_components(G_optimized)]
    cur_modularity = modularity(G_optimized, cur_components, weight=domino_weight)

    if cur_modularity >= modularity_score_objective:
        return True, best_modularity, G_optimized
    elif len(cur_components) == 0:
        return True, best_modularity, G_optimized
    else:
        optimized_connected_components = girvan_newman(G_optimized)
        cur_components = sorted(next(optimized_connected_components))
        cur_modularity = modularity(G_optimized, cur_components, weight=domino_weight)
        if cur_modularity <= best_modularity + improvement_delta:
            return True, best_modularity
        else:
            optimal_components = cur_components

            edges_to_remove = []
            for cur_edge in G_optimized.edges:
                included = False
                for n_nodes in optimal_components:
                    if cur_edge[0] in n_nodes and cur_edge[1] in n_nodes:
                        included = True
                if not included:
                    edges_to_remove.append(cur_edge)

            G_optimized.remove_edges_from(edges_to_remove)

            return False, cur_modularity, G_optimized


###########################################################################################################################

"""
Implement linear threshold models
"""

#!/usr/bin/env python
#    Copyright (C) 2004-2010 by
#    Hung-Hsuan Chen <hhchen@psu.edu>
#    All rights reserved.
#    BSD license.
#    NetworkX:http://networkx.lanl.gov/.
__author__ = """Hung-Hsuan Chen (hhchen@psu.edu)"""

__all__ = ["linear_threshold"]


def linear_threshold(G, seeds, steps=0):
    """Return the active nodes of each diffusion step by linear threshold model
    Parameters
    ----------
    G : networkx graph
        The number of nodes.
    seeds: list of nodes
        The seed nodes of the graph
    steps: int
        The number of steps to diffuse
        When steps <= 0, the model diffuses until no more nodes
        can be activated
    Return
    ------
    layer_i_nodes : list of list of activated nodes
        layer_i_nodes[0]: the seeds
        layer_i_nodes[k]: the nodes activated at the kth diffusion step
    Notes
    -----
    1. Each node is supposed to have an attribute "threshold".  If not, the
        default value is given (0.5).
    2. Each edge is supposed to have an attribute "influence".  If not, the
        default value is given (1/in_degree)
    References
    ----------
    [1] GranovetterMark. Threshold models of collective behavior.
        The American journal of sociology, 1978.
    Examples
    --------
    >>> DG = nx.DiGraph()
    >>> DG.add_edges_from([(1,2), (1,3), (1,5), (2,1), (3,2), (4,2), (4,3), \
    >>>   (4,6), (5,3), (5,4), (5,6), (6,4), (6,5)])
    >>> layers = networkx_addon.information_propagation.linear_threshold(DG, [1])
  """
    if type(G) == nx.MultiGraph or type(G) == nx.MultiDiGraph:
        raise Exception("linear_threshold() is not defined for graphs with multiedges.")

    # make sure the seeds are in the graph
    for s in seeds:
        if s not in G.nodes():
            raise Exception("seed", s, "is not in graph")

    # change to directed graph
    if not G.is_directed():
        DG = G.to_directed()
    else:
        DG = copy.deepcopy(G)

    # init thresholds
    for n in DG.nodes():
        if "threshold" not in DG.nodes[n]:
            DG.nodes[n]["threshold"] = 0.5
        elif DG.nodes[n]["threshold"] > 1:
            raise Exception(
                "node threshold:",
                DG.nodes[n]["threshold"],
                "cannot be larger than 1",
            )

    # init influences
    in_deg = DG.in_degree()
    for e in DG.edges():
        if "influence" not in DG[e[0]][e[1]]:
            DG[e[0]][e[1]]["influence"] = 1.0 / in_deg[e[1]]
        elif DG[e[0]][e[1]]["influence"] > 1:
            raise Exception(
                "edge influence:",
                DG[e[0]][e[1]]["influence"],
                "cannot be larger than 1",
            )

    # perform diffusion
    A = copy.deepcopy(seeds)
    if steps <= 0:
        # perform diffusion until no more nodes can be activated
        return _diffuse_all(DG, A)
    # perform diffusion for at most "steps" rounds only
    return _diffuse_k_rounds(DG, A, steps)


def _diffuse_all(G, A):
    layer_i_nodes = []
    layer_i_nodes.append([i for i in A])
    while True:
        len_old = len(A)
        A, activated_nodes_of_this_round = _diffuse_one_round(G, A)
        layer_i_nodes.append(activated_nodes_of_this_round)
        if len(A) == len_old:
            break
    return layer_i_nodes


def _diffuse_k_rounds(G, A, steps):
    layer_i_nodes = []
    layer_i_nodes.append([i for i in A])
    while steps > 0 and len(A) < len(G):
        len_old = len(A)
        A, activated_nodes_of_this_round = _diffuse_one_round(G, A)
        layer_i_nodes.append(activated_nodes_of_this_round)
        if len(A) == len_old:
            break
        steps -= 1
    return layer_i_nodes


def _diffuse_one_round(G, A):
    activated_nodes_of_this_round = set()
    for s in A:
        nbs = G.successors(s)
        for nb in nbs:
            if nb in A:
                continue
            active_nb = list(set(G.predecessors(nb)).intersection(set(A)))
            if _influence_sum(G, active_nb, nb) >= G.nodes[nb]["threshold"]:
                activated_nodes_of_this_round.add(nb)
    A.extend(list(activated_nodes_of_this_round))
    return A, list(activated_nodes_of_this_round)


def _influence_sum(G, froms, to):
    influence_sum = 0.0
    for f in froms:
        influence_sum += G[f][to]["influence"]
    return influence_sum
