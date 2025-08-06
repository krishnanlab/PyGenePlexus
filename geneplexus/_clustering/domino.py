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
    
    # This section is what Hao,s dta_io class does
    # This is what domino does in main lines 326 to 344
    # print(df_edge[df_edge["Node2"]=="998"])
    # make networkx grapgh object
    G = nx.from_pandas_edgelist(df_edge, source="Node1", target="Node2", edge_attr=True)
    nx.set_node_attributes(G, 0, "Score")
    nx.set_node_attributes(G, False, "Active")
    # print(G.nodes["2"])

    # make the initial slices
    if clust_weighted == True:
        domino_weight = "Weight"
    else:
        domino_weight = None
    slice_genes = nx.community.louvain_communities(G, weight=domino_weight, resolution=domino_res, seed=domino_seed)

    # ### get some extra info
    # network_genes = np.unique(np.union1d(df_edge["Node1"].to_numpy(),df_edge["Node2"].to_numpy()))
    # unique_slice_genes = list({item for sublist in slice_genes for item in sublist})
    # len_each_slice = [len(genes) for genes in slice_genes]
    # num_genes_lost = len(network_genes) - len(unique_slice_genes)
    # per_genes_lost = 100 - ((len(unique_slice_genes) / len(network_genes)) * 100)
    # print(len(network_genes))
    # print(len(unique_slice_genes))
    # print(f"The number of clusters added is {len(slice_genes)}. ")
    # print(f"The number(%) of genes lost to clustering is {num_genes_lost} ({per_genes_lost:.2f}%)")
    # print(f"The number of genes in each slice is {len_each_slice}")

    # print(input_genes)
    G = add_scores_to_G(G, input_genes)
    G_modularity = make_modularity_graph(G, slice_genes)
    G_modularity, relevant_slices, qvals = retain_relevant_slices(
        G_modularity, G, domino_slice_thresh
    )  # G_modularity and qvals maybe never called again
    ############ need to fix part of this function maybe for while loop ###############
    putative_modules = prune_slice(relevant_slices, G, domino_n_steps, clust_min_size)
    final_modules = get_final_modules(putative_modules, G, domino_module_threshold, clust_min_size)

    return final_modules


def add_scores_to_G(G, input_genes):
    # add scores and activeness to nodes if in the
    ##### the Score is never really used again
    for agene in input_genes:
        G.nodes[agene]["Score"] = 1
        G.nodes[agene]["Active"] = True
    # print(G.nodes[input_genes[0]])
    return G


def make_modularity_graph(G, slice_genes):
    # In both codes this is from prune_network_by_modularity
    # make a new full graph object for just the slice genes
    G_modules = [G.subgraph(genes) for genes in slice_genes]
    G_modularity = nx.algorithms.operators.union_all(G_modules)
    # print(G_modularity.nodes[input_genes[0]])
    return G_modularity


def retain_relevant_slices(G_modularity, G, domino_slice_thresh):
    # ##################################################################################################################
    # this is from retain relevant slices functions
    perturbed_nodes = [anode for anode in G_modularity.nodes if G_modularity.nodes[anode]["Active"]]
    ccs = [G_modularity.subgraph(c) for c in nx.connected_components(G_modularity)]
    # print(perturbed_nodes)
    # print(len(ccs))
    # set things up for the loop
    n_G_original = len(G)
    n_perturbed_nodes = len(perturbed_nodes)
    n_perturbed_nodes_in_ccs = []
    perturbed_nodes_in_ccs = []
    # print(n_G_original)
    #### n_perturbed_nodes_in_ccs can have 0 for the number of perturbed nodes in a given cc
    for i_cur_cc, cur_cc in enumerate(ccs):
        n_perturbed_nodes_in_ccs.append(
            len([cur_node for cur_node in cur_cc if G_modularity.nodes[cur_node]["Active"]]),
        )
        perturbed_nodes_in_ccs.append([cur_node for cur_node in cur_cc if G_modularity.nodes[cur_node]["Active"]])
    # print(n_perturbed_nodes_in_ccs)
    # print(perturbed_nodes_in_ccs)

    # dynamic threhold for filtering out cc with too few pertubed nodes
    # since n_G_original is large for genome wide network, the perturbation factor
    # is n of pertubed nodes in cc / n of all nodes
    perturbation_factor = min(
        0.7,
        (float(n_perturbed_nodes) / n_G_original) * (1 + 100 / n_G_original**0.5),
    )
    # when getting relevant slices, why did they not use number of pertubed nodes in cc instead of all pertubed nodes? (it is enerated in loop)
    # maybe add this but be careful between cc and ccs
    params = []
    for i_cur_cc, cur_cc in enumerate(ccs):
        params.append([cur_cc, i_cur_cc, perturbed_nodes_in_ccs[i_cur_cc], n_perturbed_nodes_in_ccs[i_cur_cc]])
    # print(params)

    ##### starting pf_filter step
    # there is a lot of switching between n_perturbed nodes and len(pertubed_nodes_in_cc)
    ######### really need to look at this section
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
            # difference icepop and dominio is icepops adds the 0.5
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
    # print(res)
    ### finishing pf_filter step
    ### back in retain_relevant_slices
    if len(res) == 0:
        G_modularity = nx.Graph()
        relevant_slices = []
        qvals = []
    else:
        large_modules, sig_scores = zip(*res)
        # print(sig_scores)
        fdr_bh_results = multipletests(
            sig_scores, method="fdr_bh", is_sorted=False, alpha=domino_slice_thresh
        )  # hao took this fdr step out
        passed_modules = [cur_cc for cur_cc, is_passed_th in zip(large_modules, fdr_bh_results[0]) if is_passed_th]
        # print(passed_modules)
        if len(passed_modules) > 0:
            G_modularity = nx.union_all(passed_modules)
            relevant_slices = [list(m.nodes) for m in passed_modules]
            qvals = [qval for qval, is_passed_th in zip(fdr_bh_results[1], fdr_bh_results[0]) if is_passed_th]
        else:
            G_modularity = nx.Graph()
            relevant_slices = []
            qvals = []
    # ice pop doesn't reurn the q or p vlaues and it looks like they might not be used again
    # I added in subset qvals to passed mdules
    # print(G_modularity)
    # print(relevant_slices)
    # print(qvals)
    # print(len(relevant_slices),len(qvals))
    # end retain_relevant_slices
    return G_modularity, relevant_slices, qvals


def prune_slice(relevant_slices, G, domino_n_steps, clust_min_size):
    ######################################################################################
    # do the prune_slice (analyze_slice in Domino)
    pruned_slices = []
    for i_cc, cc in enumerate(relevant_slices):
        G_cc = nx.subgraph(G, cc)
        nodes = list(G_cc.nodes)
        labels = {n: G_cc.nodes[n] for n in nodes}
        n_perturbed_nodes = sum([G.nodes[cur_node]["Active"] for cur_node in G.nodes])
        prize_factor = max(0, 1 - 3 * n_perturbed_nodes / float(len(G.nodes)))
        edges, edges_grid = domino_utls.run_pcst(G_cc, i_cc, labels, domino_n_steps, nodes, prize_factor)
        G_subslice = nx.Graph()
        G_subslice.add_edges_from(
            [(nodes[edges_grid[e][0]], nodes[edges_grid[e][1]]) for e in edges],
        )
        nx.set_node_attributes(G_subslice, {n: labels[n] for n in G_subslice.nodes})
        # this value range from 0.236 for np.log(10) / np.log(16000) for STRING network
        # this number do floats around modularity with some community structure
        modularity_score_objective = (
            np.log(len(G_subslice.nodes)) / np.log(len(G.nodes)) if len(G_subslice.nodes) > clust_min_size else -1
        )
        ############ need to fix part of this function maybe for while loop ###############
        subslice_after_ng, putative_modules_of_slice = domino_utls.get_putative_modules(
            G_subslice,
            improvement_delta=10**-2,
            modularity_score_objective=modularity_score_objective,
            n_cc=len(relevant_slices),
        )
        pruned_slices.append(putative_modules_of_slice)
    # flatten list of lists (changed from using reduce)
    putative_modules = [item for sublist in pruned_slices for item in sublist]
    return putative_modules


def get_final_modules(putative_modules, G, domino_module_threshold, clust_min_size):
    #####################################################################################
    # do get final moudles part
    # icepop is doing Bh and Domino does bonferroni
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
    # print(final_modules)
    return final_modules
