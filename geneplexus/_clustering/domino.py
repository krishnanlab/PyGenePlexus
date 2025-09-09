import os.path as osp

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

import geneplexus
from . import domino_utls


def domino_main(
    df_edge,
    input_genes,
    clust_min_size,
    clust_weighted,
    **kwargs,
):
    domino_res = kwargs["domino_res"]
    domino_seed = kwargs["domino_seed"]
    domino_slice_thresh = kwargs["domino_slice_thresh"]
    domino_n_steps = kwargs["domino_n_steps"]
    domino_module_threshold = kwargs["domino_module_threshold"]

    # make networkx grapgh object
    G = nx.from_pandas_edgelist(df_edge, source="Node1", target="Node2", edge_attr=True)
    nx.set_node_attributes(G, False, "Active")

    # make the initial slices
    if clust_weighted == True:
        domino_weight = "Weight"
    else:
        domino_weight = None
    slice_genes = nx.community.louvain_communities(G, weight=domino_weight, resolution=domino_res, seed=domino_seed)
    # run through rest of pipeline
    G = add_scores_to_G(G, input_genes)
    G_modularity = make_modularity_graph(G, slice_genes)
    relevant_slices = retain_relevant_slices(
        G_modularity,
        G,
        domino_slice_thresh,
    )
    putative_modules = prune_slice(relevant_slices, G, domino_n_steps, clust_min_size, domino_weight)
    final_modules = get_final_modules(putative_modules, G, domino_module_threshold, clust_min_size)

    return final_modules


def add_scores_to_G(G, input_genes):
    # input_genes are considered active genes
    for agene in input_genes:
        G.nodes[agene]["Active"] = True
    return G


def make_modularity_graph(G, slice_genes):
    # this creates a copy so G will be updated when these change
    G_modules = [G.subgraph(genes) for genes in slice_genes]
    # union_all seems to make a copy so these then won't update G
    G_modularity = nx.algorithms.operators.union_all(G_modules)
    return G_modularity


def retain_relevant_slices(G_modularity, G, domino_slice_thresh):
    perturbed_nodes = [anode for anode in G_modularity.nodes if G_modularity.nodes[anode]["Active"]]
    ccs = [G_modularity.subgraph(c) for c in nx.connected_components(G_modularity)]
    # set up loop inputs
    n_G_original = len(G)
    n_perturbed_nodes = len(perturbed_nodes)
    n_perturbed_nodes_in_ccs = []
    perturbed_nodes_in_ccs = []
    # n_perturbed_nodes_in_ccs can have 0 for the number of perturbed nodes in a given cc
    for i_cur_cc, cur_cc in enumerate(ccs):
        n_perturbed_nodes_in_ccs.append(
            len([cur_node for cur_node in cur_cc if G_modularity.nodes[cur_node]["Active"]]),
        )
        perturbed_nodes_in_ccs.append([cur_node for cur_node in cur_cc if G_modularity.nodes[cur_node]["Active"]])
    # dynamic threhold for filtering out cc with too few pertubed nodes
    # since n_G_original is large for genome wide network, the perturbation factor
    # is n of pertubed nodes in cc / n of all nodes
    perturbation_factor = min(
        0.7,
        (float(n_perturbed_nodes) / n_G_original) * (1 + 100 / n_G_original**0.5),
    )
    # maybe add this but be careful between cc and ccs
    params = []
    for i_cur_cc, cur_cc in enumerate(ccs):
        params.append([cur_cc, i_cur_cc, perturbed_nodes_in_ccs[i_cur_cc], n_perturbed_nodes_in_ccs[i_cur_cc]])
    # starting pf_filter step
    res = []
    for aparam in params:
        cur_cc, i_cur_cc, perturbed_nodes_in_cc, n_perturbed_nodes_in_cc = aparam
        if (
            len(cur_cc) < 4
            or n_perturbed_nodes == 0
            or not (
                n_perturbed_nodes_in_cc / float(len(cur_cc)) >= perturbation_factor
                or n_perturbed_nodes_in_cc / float(n_perturbed_nodes) >= 0.1
            )
        ):
            pass
        else:
            score = hypergeom.sf(
                n_perturbed_nodes_in_cc,
                n_G_original,
                n_perturbed_nodes,
                len(cur_cc),
            ) + 0.5 * hypergeom.pmf(
                n_perturbed_nodes_in_cc,
                n_G_original,
                n_perturbed_nodes,
                len(cur_cc),
            )
            res.append((cur_cc, score))
    # finishing pf_filter step
    if len(res) == 0:
        relevant_slices = []
    else:
        large_modules, sig_scores = zip(*res)
        fdr_bh_results = multipletests(
            sig_scores,
            method="fdr_bh",
            is_sorted=False,
            alpha=domino_slice_thresh,
        )  # hao took this fdr step out
        passed_modules = [cur_cc for cur_cc, is_passed_th in zip(large_modules, fdr_bh_results[0]) if is_passed_th]
        if len(passed_modules) > 0:
            relevant_slices = [list(m.nodes) for m in passed_modules]
        else:
            relevant_slices = []
    return relevant_slices


def prune_slice(relevant_slices, G, domino_n_steps, clust_min_size, domino_weight):
    pruned_slices = []
    for i_cc, cc in enumerate(relevant_slices):
        G_cc = nx.subgraph(G, cc)
        nodes = list(G_cc.nodes)
        labels = {n: G_cc.nodes[n] for n in nodes}
        n_perturbed_nodes = sum([G.nodes[cur_node]["Active"] for cur_node in G.nodes])
        prize_factor = max(0, 1 - 3 * n_perturbed_nodes / float(len(G.nodes)))
        # run_pcst subsets G_cc to only be edges important in pcst (which includes diffusion)
        edges, edges_grid = domino_utls.run_pcst(G_cc, domino_n_steps, nodes, prize_factor)
        G_subslice = nx.Graph()
        G_subslice.add_edges_from(
            [(nodes[edges_grid[e][0]], nodes[edges_grid[e][1]]) for e in edges],
        )
        nx.set_node_attributes(G_subslice, {n: labels[n] for n in G_subslice.nodes})
        # find modularity score objective
        modularity_score_objective = (
            np.log(len(G_subslice.nodes)) / np.log(len(G.nodes)) if len(G_subslice.nodes) > clust_min_size else -1
        )
        # find the putative modules
        putative_modules_of_slice = domino_utls.get_putative_modules(
            G_subslice,
            domino_weight,
            improvement_delta=10**-2,
            modularity_score_objective=modularity_score_objective,
        )
        pruned_slices.append(putative_modules_of_slice)
    # flatten list of lists
    putative_modules = [item for sublist in pruned_slices for item in sublist]
    return putative_modules


def get_final_modules(putative_modules, G, domino_module_threshold, clust_min_size):
    if len(putative_modules) == 0:
        final_modules = []
    else:
        pertubed_nodes = [cur_node for cur_node in G.nodes if G.nodes[cur_node]["Active"]]
        sig_scores = []
        for i_cur_module, cur_G_module in enumerate(putative_modules):
            pertubed_nodes_in_cc = [cur_node for cur_node in cur_G_module.nodes if G.nodes[cur_node]["Active"]]
            sig_score = hypergeom.sf(
                len(pertubed_nodes_in_cc),
                len(G.nodes),
                len(pertubed_nodes),
                len(cur_G_module.nodes),
            ) + 0.5 * hypergeom.pmf(
                len(pertubed_nodes_in_cc),
                len(G.nodes),
                len(pertubed_nodes),
                len(cur_G_module.nodes),
            )
            sig_scores.append(sig_score)
        _, adjusted_p_values, _, _ = multipletests(sig_scores, method="fdr_bh")
        # filter final modules
        module_sigs = []
        for cur_G_module, adjusted_p in zip(putative_modules, adjusted_p_values):
            if (adjusted_p <= domino_module_threshold) and (len(cur_G_module.nodes) > clust_min_size):
                module_sigs.append((cur_G_module, adjusted_p))
        module_sigs = sorted(module_sigs, key=lambda a: a[1])
        final_modules = [list(a[0].nodes) for a in module_sigs]
    return final_modules
