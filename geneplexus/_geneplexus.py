import os.path as osp
import pathlib
import shutil
from typing import Any
from typing import Dict
from typing import Optional

import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine
from scipy.stats import hypergeom
from scipy.stats import norm
from scipy.stats import rankdata
from scipy.stats import zscore
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

from . import util
from ._clustering.domino import domino_main
from ._clustering.louvain import louvain_main
from ._config import logger
from ._config.config import DEFAULT_LOGREG_KWARGS


def _initial_id_convert(input_genes, file_loc, species):
    # load all the possible conversion dictionaries
    convert_types = ["ENSG", "Symbol", "ENSP", "ENST"]
    all_convert_dict = {}
    for anIDtype in convert_types:
        convert_tmp = util.load_geneid_conversion(file_loc, species, anIDtype, "Entrez", upper=True)
        all_convert_dict[anIDtype] = convert_tmp
    # find all Entrez genes in networks and conversion files
    all_entrez_genes = np.array([])
    for anet in util.get_all_net_types(file_loc, species):
        if (species == "Zebrafish") and (anet == "BioGRID"):
            continue
        entrez_genes_tmp = util.load_node_order(file_loc, species, anet)
        all_entrez_genes = np.union1d(all_entrez_genes, entrez_genes_tmp)
    entrez_genes_tmp = util.load_geneid_conversion(file_loc, species, "Entrez", "Symbol")
    all_entrez_genes = np.union1d(all_entrez_genes, np.array(list(entrez_genes_tmp)))
    entrez_to_name = util.load_geneid_conversion(file_loc, species, "Entrez", "Name")
    all_entrez_genes = np.union1d(all_entrez_genes, np.array(list(entrez_to_name)))

    # make some place holder arrays
    convert_ids = []  # This will be a flat list for Entrez IDs to use as positives
    convert_out = []  # This will be a list of lists that will be used to tell user the conversions made
    for agene in input_genes:
        try:
            agene_int = int(agene)
            if agene in all_entrez_genes:
                agene_name = ", ".join(entrez_to_name[agene])
                convert_out.append([agene, agene, agene_name])
                convert_ids.append(agene)
            else:
                convert_out.append(
                    [agene, f"Not in Our List of {species} Entrez Genes", None],
                )
        except ValueError:
            converted_gene: Optional[str] = None
            converted_gene_name: Optional[str] = None
            for anIDtype in convert_types:
                if agene in all_convert_dict[anIDtype]:
                    convert_ids.extend(all_convert_dict[anIDtype][agene])
                    converted_gene = ", ".join(all_convert_dict[anIDtype][agene])
                    all_to_name = []
                    for ent_gene in all_convert_dict[anIDtype][agene]:
                        all_to_name = all_to_name + entrez_to_name[ent_gene]
                    converted_gene_name = ", ".join(all_to_name)
                    logger.debug(f"Found mapping ({anIDtype}) {agene} -> {all_convert_dict[anIDtype][agene]}")
                    break
            convert_out.append(
                [
                    agene,
                    converted_gene or "Could Not be mapped to Entrez",
                    converted_gene_name or None,
                ],
            )

    column_names = ["Original ID", "Entrez ID", "Gene Name"]
    df_convert_out = pd.DataFrame(convert_out, columns=column_names).astype(str)

    return convert_ids, df_convert_out


def _make_validation_df(df_convert_out, file_loc, species):
    table_summary = []
    input_count = df_convert_out.shape[0]
    converted_genes = df_convert_out["Entrez ID"].to_numpy()
    for anet in util.get_all_net_types(file_loc, species):
        if (species == "Zebrafish") and (anet == "BioGRID"):
            continue
        net_genes = util.load_node_order(file_loc, species, anet)
        df_tmp = df_convert_out[df_convert_out["Entrez ID"].isin(net_genes)]
        pos_genes_in_net = np.intersect1d(converted_genes, net_genes)
        table_row = {"Network": anet, "NetworkGenes": len(net_genes), "PositiveGenes": len(pos_genes_in_net)}
        table_summary.append(dict(table_row))
        tmp_ins = np.full(len(converted_genes), "N", dtype=str)
        tmp_ins[df_tmp.index.to_numpy()] = "Y"
        df_convert_out[f"In {anet}?"] = tmp_ins

    return df_convert_out, table_summary, input_count


def _generate_clusters(
    file_loc,
    species,
    net_type,
    input_genes,
    clust_method,
    clust_min_size,
    clust_weighted,
    **kwargs,
):
    logger.info("Generating Clusters")
    # Load network as edge list dataframe
    filepath = osp.join(file_loc, f"Edgelist__{species}__{net_type}.edg")
    if net_type == "BioGRID":
        df_edge = pd.read_csv(filepath, sep="\t", header=None, names=["Node1", "Node2"])
        df_edge["Weight"] = 1.0
    else:
        df_edge = pd.read_csv(filepath, sep="\t", header=None, names=["Node1", "Node2", "Weight"])
        df_edge["Weight"] = df_edge["Weight"].astype(float)
    df_edge = df_edge.astype({"Node1": str, "Node2": str})
    # edgelist used for clustering can be subset of genes in data and node_order
    # make sure those input geens are removed
    df_edge_genes = np.unique(np.union1d(df_edge["Node1"].to_numpy(), df_edge["Node2"].to_numpy()))
    bad_input_genes = np.setdiff1d(input_genes, df_edge_genes).tolist()
    input_genes = np.intersect1d(input_genes, df_edge_genes).tolist()
    logger.info(
        f"Input genes removed that weren't in network being clustered (thresholded version of the full network) {bad_input_genes}",
    )
    # iteratively run through clustering
    if clust_method == "louvain":
        final_clusters = louvain_main(
            df_edge,
            input_genes,
            clust_min_size,
            clust_weighted,
            **kwargs,
        )
    elif clust_method == "domino":
        final_clusters = domino_main(
            df_edge,
            input_genes,
            clust_min_size,
            clust_weighted,
            **kwargs,
        )
    return final_clusters


def _get_genes_in_network(file_loc, species, net_type, convert_ids):
    net_genes = util.load_node_order(file_loc, species, net_type)
    convert_ids = np.array(convert_ids, dtype=str)
    pos_genes_in_net = np.intersect1d(convert_ids, net_genes)
    genes_not_in_net = np.setdiff1d(convert_ids, net_genes)
    return pos_genes_in_net, genes_not_in_net, net_genes


def _get_negatives(file_loc, species, net_type, gsc, pos_genes_in_net, user_negatives):
    gsc_full = util.load_gsc(file_loc, species, gsc, net_type)
    uni_genes = np.array(gsc_full["Universe"])
    gsc_terms = gsc_full["Term_Order"]
    M = len(uni_genes)
    N = len(pos_genes_in_net)
    neutral_gene_info = {}
    genes_to_remove = pos_genes_in_net
    for akey in gsc_terms:
        if user_negatives == None:
            n_set = gsc_full[akey]["Genes"]
        else:
            n_set = np.setdiff1d(gsc_full[akey]["Genes"], user_negatives).tolist()
        n = len(n_set)
        k = len(np.intersect1d(pos_genes_in_net, n_set))
        pval = hypergeom.sf(k - 1, M, n, N)
        if pval < 0.05:
            genes_to_remove = np.union1d(genes_to_remove, n_set)
            neutral_gene_info[akey] = {
                "Name": gsc_full[akey]["Name"],
                "Task": gsc_full[akey]["Task"],
                "Genes": n_set,
            }
    neutral_gene_info["All Neutrals"] = np.setdiff1d(genes_to_remove, pos_genes_in_net).tolist()
    negative_genes = np.setdiff1d(uni_genes, genes_to_remove)
    return negative_genes, neutral_gene_info


def _run_sl(
    file_loc,
    sp_trn,
    net_type,
    features,
    pos_genes_in_net,
    negative_genes,
    net_genes,
    logreg_kwargs,
    min_num_pos_cv,
    num_folds,
    null_val,
    random_state,
    cross_validate,
    scale,
):
    logreg_kwargs_defaults = DEFAULT_LOGREG_KWARGS
    if isinstance(logreg_kwargs, dict):
        logreg_kwargs_defaults.update(logreg_kwargs)
    else:
        logger.warning(f"logreg_kwargs not a dictionary or None, using defaults")
    # train the model
    pos_inds = [np.where(net_genes == agene)[0][0] for agene in pos_genes_in_net]
    neg_inds = [np.where(net_genes == agene)[0][0] for agene in negative_genes]
    data = util.load_gene_features(file_loc, sp_trn, features, net_type)
    if scale:
        std_scale = StandardScaler().fit(data)
        data = std_scale.transform(data)
    Xdata = data[pos_inds + neg_inds, :]
    ydata = np.array([1] * len(pos_inds) + [0] * len(neg_inds))
    clf = LogisticRegression(**logreg_kwargs_defaults)
    clf.fit(Xdata, ydata)
    mdl_weights = np.squeeze(clf.coef_)
    # validate the model
    avgps = [null_val] * num_folds
    if not cross_validate:
        logger.info("Skipping cross validation.")
    elif len(pos_genes_in_net) < min_num_pos_cv:
        logger.warning(
            "Insufficient number of positive genes for cross validation: "
            f"{len(pos_genes_in_net)} ({min_num_pos_cv} needed). Skipping cross "
            f"validation and filling values with {null_val}",
        )
    else:
        logger.info("Performing cross validation.")
        avgps = []
        skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=random_state)
        for trn_inds, tst_inds in skf.split(Xdata, ydata):
            clf_cv = LogisticRegression(**logreg_kwargs_defaults)
            clf_cv.fit(Xdata[trn_inds], ydata[trn_inds])
            probs_cv = clf_cv.predict_proba(Xdata[tst_inds])[:, 1]
            avgp = average_precision_score(ydata[tst_inds], probs_cv)
            num_tst_pos = np.sum(ydata[tst_inds])
            prior = num_tst_pos / Xdata[tst_inds].shape[0]
            log2_prior = float(np.log2(avgp / prior))
            avgps.append(log2_prior)
        logger.info(f"{avgps=}")
        logger.info(f"{np.median(avgps)=:.2f}")
        logger.info(f"{np.mean(avgps)=:.2f}")
    if scale:
        return mdl_weights, avgps, clf, std_scale
    else:
        return mdl_weights, avgps, clf, None


def _get_predictions(
    file_loc,
    sp_res,
    features,
    net_type,
    scale,
    std_scale,
    clf,
):
    # do predictions in the target species
    data = util.load_gene_features(file_loc, sp_res, features, net_type)
    if scale:
        data = std_scale.transform(data)
    probs = clf.predict_proba(data)[:, 1]
    return probs


def _make_prob_df(file_loc, sp_trn, sp_res, net_type, probs, pos_genes_in_net, negative_genes):
    Entrez_to_Symbol = util.load_geneid_conversion(file_loc, sp_res, "Entrez", "Symbol")
    Entrez_to_Name = util.load_geneid_conversion(file_loc, sp_res, "Entrez", "Name")
    net_genes = util.load_node_order(file_loc, sp_res, net_type)
    if sp_trn != sp_res:
        biomart_orthos = util.load_biomart(file_loc, sp_trn, sp_res)
        pos_genes_tmp = [biomart_orthos[item] for item in pos_genes_in_net if item in biomart_orthos]
        neg_genes_tmp = [biomart_orthos[item] for item in negative_genes if item in biomart_orthos]
    else:
        pos_genes_tmp = pos_genes_in_net
        neg_genes_tmp = negative_genes
    prob_results = []
    for idx in range(len(net_genes)):
        if net_genes[idx] in pos_genes_tmp:
            class_label = "P"
            novel_label = "Known"
        elif net_genes[idx] in neg_genes_tmp:
            class_label = "N"
            novel_label = "Novel"
        else:
            class_label = "U"
            novel_label = "Novel"
        syms_tmp = util.mapgene(net_genes[idx], Entrez_to_Symbol)
        name_tmp = util.mapgene(net_genes[idx], Entrez_to_Name)
        prob_results.append([net_genes[idx], syms_tmp, name_tmp, novel_label, class_label, probs[idx]])
    df_col_names = ["Entrez", "Symbol", "Name", "Known/Novel", "Class-Label", "Probability"]
    df_probs = pd.DataFrame(prob_results, columns=df_col_names)
    df_probs = df_probs.astype({"Entrez": str, "Probability": float})
    df_probs = df_probs.sort_values(by=["Probability"], ascending=False).reset_index(drop=True)
    z = zscore(df_probs["Probability"].to_numpy())
    p = norm.sf(z)
    rejects, padjusts, b, c = multipletests(p, method="bonferroni", is_sorted=True)
    df_probs["Z-score"] = z
    df_probs["P-adjusted"] = padjusts
    df_probs["Rank"] = rankdata(1 / (df_probs["Probability"].to_numpy() + 1e-9), method="min")
    return df_probs


def _make_sim_dfs(file_loc, mdl_weights, species, gsc, net_type, features):
    # make task conversion
    task_convert = {"Mondo": "Disease", "Monarch": "Phenotype", "GO": "Biological Process"}
    weights_dict = util.load_pretrained_weights(file_loc, species, gsc, net_type, features)
    gsc_full = util.load_gsc(file_loc, species, gsc, net_type)
    gsc_terms = gsc_full["Term_Order"]
    mdl_sims = []
    for idx, aset in enumerate(gsc_terms):
        cos_sim = 1 - cosine(weights_dict[aset]["Weights"], mdl_weights)
        mdl_sims.append(cos_sim)
    mdl_sims = np.array(mdl_sims)
    z = zscore(mdl_sims)
    results_tmp = []
    for idx2, termID_tmp in enumerate(gsc_terms):
        if gsc_full[termID_tmp]["Task"] in task_convert:
            Task = task_convert[gsc_full[termID_tmp]["Task"]]
        else:
            Task = gsc_full[termID_tmp]["Task"]
        ID_tmp = termID_tmp
        Name_tmp = weights_dict[termID_tmp]["Name"]
        mdl_sim_tmp = mdl_sims[idx2]
        z_tmp = z[idx2]
        results_tmp.append([Task, ID_tmp, Name_tmp, mdl_sim_tmp, z_tmp])
    df_sim = pd.DataFrame(results_tmp, columns=["Task", "ID", "Name", "Similarity", "Z-score"]).sort_values(
        by=["Similarity"],
        ascending=False,
    )
    p = norm.sf(df_sim["Z-score"].to_numpy())
    rejects, padjusts, b, c = multipletests(p, method="bonferroni", is_sorted=True)
    df_sim["P-adjusted"] = padjusts
    df_sim["Rank"] = rankdata(-1 * (df_sim["Similarity"].to_numpy() + 1e-9), method="min")
    return df_sim, weights_dict


def _make_small_edgelist(file_loc, df_probs, species, net_type, num_nodes=50):
    # This will set the max number of genes to look at to a given number
    # Load network as edge list dataframe
    filepath = osp.join(file_loc, f"Edgelist__{species}__{net_type}.edg")
    if net_type == "BioGRID":
        df_edge = pd.read_csv(filepath, sep="\t", header=None, names=["Node1", "Node2"])
    else:
        df_edge = pd.read_csv(filepath, sep="\t", header=None, names=["Node1", "Node2", "Weight"])
    df_edge = df_edge.astype({"Node1": str, "Node2": str})
    # Take subgraph induced by top genes
    top_genes = df_probs["Entrez"].to_numpy()[:num_nodes]
    df_edge = df_edge[(df_edge["Node1"].isin(top_genes)) & (df_edge["Node2"].isin(top_genes))]
    if net_type == "BioGRID":
        df_edge["Weight"] = [1.0] * df_edge.shape[0]
    df_edge = df_edge.sort_values(by=["Weight", "Node1", "Node2"], ascending=[False, True, True])
    genes_in_edge = np.union1d(df_edge["Node1"].unique(), df_edge["Node2"].unique())
    isolated_genes = np.setdiff1d(top_genes, genes_in_edge).tolist()
    isolated_genes = [str(item) for item in isolated_genes]
    # Convert to gene symbol
    Entrez_to_Symbol = util.load_geneid_conversion(file_loc, species, "Entrez", "Symbol")
    replace_dict = {gene: util.mapgene(gene, Entrez_to_Symbol) for gene in genes_in_edge}
    isolated_genes_sym = [util.mapgene(gene, Entrez_to_Symbol) for gene in isolated_genes]
    df_edge_sym = df_edge.replace(to_replace=replace_dict)
    df_edge_sym = df_edge_sym.sort_values(by=["Weight", "Node1", "Node2"], ascending=[False, True, True])
    return df_edge, isolated_genes, df_edge_sym, isolated_genes_sym


def _alter_validation_df(df_convert_out, pos_genes_for_model, net_type):
    df_convert_out_subset = df_convert_out[["Original ID", "Entrez ID", "Gene Name", f"In {net_type}?"]]
    df_convert_out_subset = df_convert_out_subset[df_convert_out_subset["Entrez ID"].isin(pos_genes_for_model)]
    return df_convert_out_subset


def _save_class(gp, output_dir, save_type, zip_output, overwrite):
    output_dir = util.suffix_dir(output_dir, overwrite=overwrite)
    util._save_results(gp, output_dir, save_type)
    # Optionally zip the result directory
    if zip_output:
        outpath = pathlib.Path(output_dir)
        zip_outpath = f"{outpath}.zip"
        logger.info("Zipping output files")
        shutil.make_archive(zip_outpath[:-4], "zip", outpath.parent, outpath.name)
        shutil.rmtree(output_dir)
        logger.info(f"Removing temporary directory {output_dir}")
        logger.info(f"Done! Results saved to {zip_outpath}")
    else:
        logger.info(f"Done! Results saved to {output_dir}")
