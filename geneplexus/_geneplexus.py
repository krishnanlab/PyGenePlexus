import os.path as osp
from typing import Optional

import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine
from scipy.stats import hypergeom
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from . import config
from . import util
from ._config import logger


def initial_ID_convert(input_genes, file_loc):
    # load all the possible conversion dictionaries
    convert_types = ["ENSG", "Symbol", "ENSP", "ENST"]
    all_convert_dict = {}
    for anIDtype in convert_types:
        convert_tmp = util.load_geneid_conversion(file_loc, anIDtype, "Entrez", upper=True)
        all_convert_dict[anIDtype] = convert_tmp

    # make some place holder arrays
    convert_IDs = []  # This will be a flat list for Entrez IDs to use as positives
    convert_out = []  # This will be a list of lists that will be used to tell user the conversions made
    for agene in input_genes:
        try:
            agene_int = int(agene)
            convert_out.append([agene_int, agene_int])
            convert_IDs.append(agene_int)
        except ValueError:
            converted_gene: Optional[str] = None
            for anIDtype, conversion_map in convert_types.items():
                if agene in conversion_map:
                    convert_IDs.extend(conversion_map[agene])
                    converted_gene = ", ".join(conversion_map[agene])
                    logger.debug(f"Found mapping ({anIDtype}) {agene} -> {conversion_map[agene]}")
                    break
            convert_out.append([agene, converted_gene or "Could Not be mapped to Entrez"])

    column_names = ["Original ID", "Entrez ID"]
    df_convert_out = pd.DataFrame(convert_out, columns=column_names).astype(str)

    return convert_IDs, df_convert_out


def make_validation_df(df_convert_out, file_loc):
    table_summary = []
    input_count = df_convert_out.shape[0]
    converted_genes = df_convert_out["Entrez ID"].to_numpy()

    for anet in config.ALL_NETWORKS:
        net_genes = util.load_node_order(file_loc, anet)
        df_tmp = df_convert_out[df_convert_out["Entrez ID"].isin(net_genes)]
        pos_genes_in_net = np.intersect1d(converted_genes, net_genes)
        table_row = {"Network": anet, "NetworkGenes": len(net_genes), "PositiveGenes": len(pos_genes_in_net)}
        table_summary.append(dict(table_row))
        tmp_ins = np.full(len(converted_genes), "N", dtype=str)
        tmp_ins[df_tmp.index.to_numpy()] = "Y"
        df_convert_out[f"In {anet}?"] = tmp_ins

    return df_convert_out, table_summary, input_count


def get_genes_in_network(file_loc, net_type, convert_IDs):
    net_genes = util.load_node_order(file_loc, net_type)
    convert_IDs_array = np.array(convert_IDs, dtype=str)
    pos_genes_in_net = np.intersect1d(convert_IDs_array, net_genes)
    genes_not_in_net = np.setdiff1d(convert_IDs_array, net_genes)
    return pos_genes_in_net, genes_not_in_net, net_genes


def get_negatives(file_loc, net_type, GSC, pos_genes_in_net):
    uni_genes = util.load_genes_universe(file_loc, GSC, net_type)
    good_sets = util.load_gsc(file_loc, GSC, net_type)
    M = len(uni_genes)
    N = len(pos_genes_in_net)
    genes_to_remove = pos_genes_in_net
    for akey in good_sets:
        n = len(good_sets[akey]["Genes"])
        k = len(np.intersect1d(pos_genes_in_net, good_sets[akey]["Genes"]))
        pval = hypergeom.sf(k - 1, M, n, N)
        if pval < 0.05:
            genes_to_remove = np.union1d(genes_to_remove, good_sets[akey]["Genes"])
    negative_genes = np.setdiff1d(uni_genes, genes_to_remove)
    return negative_genes


def run_SL(file_loc, net_type, features, pos_genes_in_net, negative_genes, net_genes):
    pos_inds = [np.where(net_genes == agene)[0][0] for agene in pos_genes_in_net]
    neg_inds = [np.where(net_genes == agene)[0][0] for agene in negative_genes]
    data = util.load_gene_features(file_loc, features, net_type)
    std_scale = StandardScaler().fit(data)
    data = std_scale.transform(data)
    Xdata = data[pos_inds + neg_inds, :]
    ydata = np.array([1] * len(pos_inds) + [0] * len(neg_inds))
    clf = LogisticRegression(max_iter=10000, solver="lbfgs", penalty="l2", C=1.0)
    clf.fit(Xdata, ydata)
    mdl_weights = np.squeeze(clf.coef_)
    probs = clf.predict_proba(data)[:, 1]

    if len(pos_genes_in_net) < 15:
        logger.warning(
            "Not enough genes supplied for cross validation ({len(pos_genes_in_net)}), skipping cross validation.",
        )
        avgps = [-10, -10, -10]
    else:
        avgps = []
        n_folds = 3
        skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=None)
        for trn_inds, tst_inds in skf.split(Xdata, ydata):
            clf_cv = LogisticRegression(max_iter=10000, solver="lbfgs", penalty="l2", C=1.0)
            clf_cv.fit(Xdata[trn_inds], ydata[trn_inds])
            probs_cv = clf_cv.predict_proba(Xdata[tst_inds])[:, 1]
            avgp = average_precision_score(ydata[tst_inds], probs_cv)
            num_tst_pos = np.sum(ydata[tst_inds])
            prior = num_tst_pos / Xdata[tst_inds].shape[0]
            log2_prior = np.log2(avgp / prior)
            avgps.append(log2_prior)
        logger.info(f"{avgps=}")
        logger.info(f"{np.median(avgps)=:.2f}")
        logger.info(f"{np.mean(avgps)=:.2f}")
    return mdl_weights, probs, avgps


def make_prob_df(file_loc, net_genes, probs, pos_genes_in_net, negative_genes):
    Entrez_to_Symbol = util.load_geneid_conversion(file_loc, "Entrez", "Symbol")
    Entrez_to_Name = util.load_geneid_conversion(file_loc, "Entrez", "Name")
    prob_results = []
    for idx in range(len(net_genes)):
        if net_genes[idx] in pos_genes_in_net:
            class_label = "P"
            novel_label = "Known"
        elif net_genes[idx] in negative_genes:
            class_label = "N"
            novel_label = "Novel"
        else:
            class_label = "U"
            novel_label = "Novel"
        try:
            syms_tmp = "/".join(Entrez_to_Symbol[net_genes[idx]])  # allows for multimapping
        except KeyError:
            syms_tmp = "N/A"
        try:
            name_tmp = "/".join(Entrez_to_Name[net_genes[idx]])  # allows for multimapping
        except KeyError:
            name_tmp = "N/A"
        prob_results.append([net_genes[idx], syms_tmp, name_tmp, probs[idx], novel_label, class_label])
    df_probs = pd.DataFrame(
        prob_results,
        columns=["Entrez", "Symbol", "Name", "Probability", "Known/Novel", "Class-Label"],
    )
    df_probs = df_probs.astype({"Entrez": str, "Probability": float})
    df_probs = df_probs.sort_values(by=["Probability"], ascending=False)
    df_probs["Rank"] = rankdata(1 / (df_probs["Probability"].to_numpy() + 1e-9), method="min")
    return df_probs


def make_sim_dfs(file_loc, mdl_weights, GSC, net_type, features):
    dfs_out = []
    for target_set in ["GO", "DisGeNet"]:
        weights_dict = util.load_pretrained_weights(file_loc, target_set, net_type, features)
        if target_set == "GO":
            weights_dict_GO = weights_dict
        if target_set == "DisGeNet":
            weights_dict_Dis = weights_dict
        order = util.load_correction_order(file_loc, target_set, net_type)
        cor_mat = util.load_correction_mat(file_loc, GSC, target_set, net_type, features)
        add_row = np.zeros((1, len(order)))
        for idx, aset in enumerate(order):
            cos_sim = 1 - cosine(weights_dict[aset]["Weights"], mdl_weights)
            add_row[0, idx] = cos_sim
        cor_mat = np.concatenate((cor_mat, add_row), axis=0)
        last_row = cor_mat[-1, :]
        zq = np.maximum(0, (last_row - np.mean(last_row)) / np.std(last_row))
        zs = np.maximum(0, (last_row - np.mean(cor_mat, axis=0)) / np.std(cor_mat, axis=0))
        z = np.sqrt(zq ** 2 + zs ** 2)
        results_tmp = []
        for idx2, termID_tmp in enumerate(order):
            ID_tmp = termID_tmp
            Name_tmp = weights_dict[termID_tmp]["Name"]
            z_tmp = z[idx2]
            results_tmp.append([ID_tmp, Name_tmp, z_tmp])
        df_tmp = pd.DataFrame(results_tmp, columns=["ID", "Name", "Similarity"]).sort_values(
            by=["Similarity"],
            ascending=False,
        )
        df_tmp["Rank"] = rankdata(1 / (df_tmp["Similarity"].to_numpy() + 1e-9), method="min")
        dfs_out.append(df_tmp)
    return dfs_out[0], dfs_out[1], weights_dict_GO, weights_dict_Dis


def make_small_edgelist(file_loc, df_probs, net_type, num_nodes=50):
    # This will set the max number of genes to look at to a given number
    # Load network as edge list dataframe
    filepath = osp.join(file_loc, f"Edgelist_{net_type}.edg")
    if net_type == "BioGRID":
        df_edge = pd.read_csv(filepath, sep="\t", header=None, names=["Node1", "Node2"])
        df_edge["Weight"] = 1
    else:
        df_edge = pd.read_csv(filepath, sep="\t", header=None, names=["Node1", "Node2", "Weight"])
    df_edge = df_edge.astype({"Node1": str, "Node2": str})

    # Take subgraph induced by top genes
    top_genes = df_probs["Entrez"].to_numpy()[0:num_nodes]
    df_edge = df_edge[(df_edge["Node1"].isin(top_genes)) & (df_edge["Node2"].isin(top_genes))]
    genes_in_edge = np.union1d(df_edge["Node1"].unique(), df_edge["Node2"].unique())
    isolated_genes = np.setdiff1d(top_genes, genes_in_edge).tolist()

    # Convert to gene symbol
    Entrez_to_Symbol = util.load_geneid_conversion(file_loc, "Entrez", "Symbol")
    replace_dict = {gene: util.get_symbol(gene, Entrez_to_Symbol) for gene in genes_in_edge}
    isolated_genes_sym = [util.get_symbol(gene, Entrez_to_Symbol) for gene in isolated_genes]
    df_edge_sym = df_edge.replace(to_replace=replace_dict)

    return df_edge, isolated_genes, df_edge_sym, isolated_genes_sym


def alter_validation_df(df_convert_out, table_summary, net_type):
    df_convert_out_subset = df_convert_out[["Original ID", "Entrez ID", f"In {net_type}?"]]
    network = next((item for item in table_summary if item["Network"] == net_type), None)
    positive_genes = network.get("PositiveGenes")
    return df_convert_out_subset, positive_genes
