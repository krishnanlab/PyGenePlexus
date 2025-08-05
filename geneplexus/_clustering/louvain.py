import networkx as nx


def louvain_main(
    df_edge,
    input_genes,
    clust_min_size,
    louvain_max_size,
    louvain_max_tries,
    louvain_res,
    clust_weighted,
    louvain_seed,
):
    for num_try in range(louvain_max_tries):
        ######## logger.info(f"On clustering try {clus_try + 1}")
        if num_try == 0:
            final_clusters, large_clusters = louvain_cluster(
                df_edge,
                [input_genes],
                [],
                clust_min_size,
                louvain_max_size,
                louvain_res,
                clust_weighted,
                louvain_seed,
            )
        else:
            final_clusters, large_clusters = louvain_cluster(
                df_edge,
                large_clusters,
                final_clusters,
                clust_min_size,
                louvain_max_size,
                louvain_res,
                clust_weighted,
                louvain_seed,
            )
        if len(large_clusters) == 0:
            break
    final_clusters = final_clusters + large_clusters  # add back in large clusters if couldn't be made smaller
    return final_clusters


def louvain_cluster(
    df_edge,
    sets_to_cluster,
    final_clusters,
    clust_min_size,
    louvain_max_size,
    louvain_res,
    clust_weighted,
    louvain_seed,
):
    if clust_weighted == True:
        louvain_weight = "Weight"
    else:
        louvain_weight = None
    for aset in sets_to_cluster:
        df_edge_tmp = df_edge[(df_edge["Node1"].isin(aset)) & (df_edge["Node2"].isin(aset))]
        G = nx.from_pandas_edgelist(df_edge_tmp, source="Node1", target="Node2", edge_attr=True)
        clusters = nx.community.louvain_communities(G, weight=louvain_weight, resolution=louvain_res, seed=louvain_seed)
        large_clusters = []
        for idx, aclus in enumerate(clusters):
            if len(aclus) >= clust_min_size and len(aclus) <= louvain_max_size:
                final_clusters.append(list(aclus))
            elif len(aclus) > louvain_max_size:
                large_clusters.append(list(aclus))
    return final_clusters, large_clusters
