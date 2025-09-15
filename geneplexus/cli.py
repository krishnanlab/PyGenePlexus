"""Command line interface for the GenePlexus pipeline."""
import argparse
import atexit
import json
import os
import os.path as osp
import pathlib
import shutil

import numpy as np
import pandas as pd

from . import config
from ._config import logger
from .geneplexus import GenePlexus
from .util import format_choices
from .util import normexpand
from .util import read_gene_list

os.environ["COLUMNS"] = "100"  # for CLI help page wrap line


def parse_args() -> argparse.Namespace:
    """Parse arguments from command line."""
    parser = argparse.ArgumentParser(
        description="Run the GenePlexus pipline on a input gene list.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=argparse.SUPPRESS,
    )

    ### main set of arguements ###

    parser.add_argument(
        "-i",
        "--input_file",
        metavar="",
        required=True,
        help="Input gene list (.txt) file.",
    )

    parser.add_argument(
        "-d",
        "--gene_list_delimiter",
        default="newline",
        metavar="",
        help="Delimiter used in the gene list. Use 'newline' if the genes are "
        "separated by new line, and use 'tab' if the genes are seperate by "
        "tabs. If not newline or tab, will use argument directly, so /t, /n, ,",
    )

    parser.add_argument(
        "-fl",
        "--file_loc",
        default=config.DEFAULT_PARAMETERS["file_loc"],
        metavar="",
        help="Directory in which the data are stored, if set to None, then use "
        "the default data directory ~/.data/geneplexus",
    )

    parser.add_argument(
        "-n",
        "--network",
        default="STRING",
        metavar="",
        help=f"Network to use. {format_choices(config.ALL_NETWORKS)}",
    )

    parser.add_argument(
        "-f",
        "--feature",
        default="SixSpeciesN2V",
        metavar="",
        help=f"Types of feature to use. {format_choices(config.ALL_FEATURES)}",
    )

    parser.add_argument(
        "-s1",
        "--sp_trn",
        default="Human",
        metavar="",
        help=f"Species of training data {format_choices(config.ALL_SPECIES)}",
    )

    parser.add_argument(
        "-s2",
        "--sp_res",
        default="Mouse",
        metavar="",
        help=f"Species of results data {format_choices(config.ALL_SPECIES)}. "
        "If more than one species make comma seaprated.",
    )

    parser.add_argument(
        "-g1",
        "--gsc_trn",
        default="GO",
        metavar="",
        help=f"Geneset collection used to generate negatives. {format_choices(config.ALL_GSCS)}",
    )

    parser.add_argument(
        "-g2",
        "--gsc_res",
        default="GO",
        metavar="",
        help=f"Geneset collection used for model similarities. {format_choices(config.ALL_GSCS)}. "
        "If more than one gsc can be comma spearated.",
    )

    parser.add_argument(
        "-od",
        "--output_dir",
        default="result/",
        metavar="",
        help="Output directory with respect to the repo root directory.",
    )

    parser.add_argument(
        "-in",
        "--input_negatives",
        default=None,
        metavar="",
        help="Input negative gene list (.txt) file.",
    )

    ### pipeline control arguements ###

    parser.add_argument(
        "-ad",
        "--auto_download",
        action="store_true",
        help="When added turns on autodownloader which is off by default.",
    )

    parser.add_argument(
        "--clear-data",
        action="store_true",
        help="Clear data directory and exit.",
    )

    parser.add_argument(
        "--do_clustering",
        action="store_true",
        help="Do clustering step.",
    )

    parser.add_argument(
        "--skip-mdl-sim",
        action="store_true",
        help="Skip model similarity computation",
    )

    parser.add_argument(
        "--skip-sm-edgelist",
        action="store_true",
        help="Skip making small edgelist.",
    )

    ### Class methods arguements ###

    parser.add_argument(
        "-l",
        "--log_level",
        default="INFO",
        metavar="",
        help=f"Logging level. {format_choices(config.LOG_LEVELS)}",
    )

    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress log messages (same as setting log_level to CRITICAL).",
    )

    parser.add_argument(
        "-cmin",
        "--clust_min_size",
        default=5,
        metavar="",
        type=int,
        help="Minimum size of clusters allowed.",
    )

    parser.add_argument(
        "-cmax",
        "--clust_max_size",
        default=70,
        metavar="",
        type=int,
        help="Maximum size of clusters allowed.",
    )

    parser.add_argument(
        "-ctries",
        "--clust_max_tries",
        default=3,
        metavar="",
        type=int,
        help="Number of times to try to sub-cluster large clusters.",
    )

    parser.add_argument(
        "-cres",
        "--clust_res",
        default=1,
        metavar="",
        type=int,
        help="Cluster resolution parameter.",
    )

    parser.add_argument(
        "-cweight",
        "--clust_unweighted",
        action="store_false",
        help="If set, will not use cluster weights.",
    )

    parser.add_argument(
        "-flk",
        "--fit_logreg_kwargs",
        default=None,
        metavar="",
        type=json.loads,
        help="Logistic regression leyword arguments.",
    )

    parser.add_argument(
        "-fs",
        "--fit_scale",
        action="store_true",
        help="If set, will scale input data. See docs for more info of when this is good to do.",
    )

    parser.add_argument(
        "-fmnp",
        "--fit_min_num_pos",
        default=5,
        metavar="",
        type=int,
        help="Number of genes needed to fit a model.",
    )

    parser.add_argument(
        "-fmnpcv",
        "--fit_min_num_pos_cv",
        default=15,
        metavar="",
        type=int,
        help="Number of genes needed to do cross validation.",
    )

    parser.add_argument(
        "-fnf",
        "--fit_num_folds",
        default=3,
        metavar="",
        type=int,
        help="Number of genes needed to do cross validation.",
    )

    parser.add_argument(
        "-fnv",
        "--fit_null_val",
        default=None,
        metavar="",
        type=float,
        help="Value to use when CV can't be done.",
    )

    parser.add_argument(
        "-frs",
        "--fit_random_state",
        default=0,
        metavar="",
        type=int,
        help="Random state value to use when fitting.",
    )

    parser.add_argument(
        "-fscv",
        "--fit_skip_cross_validate",
        action="store_false",
        help="If set, will not try to do CV.",
    )

    parser.add_argument(
        "-s",
        "--small_edgelist_num_nodes",
        default=50,
        metavar="",
        type=int,
        help="Number of nodes in the small edgelist.",
    )

    parser.add_argument(
        "-st",
        "--save_type",
        default="all",
        metavar="",
        type=str,
        help="Which files to save.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing result directory if set.",
    )

    parser.add_argument(
        "-z",
        "--zip-output",
        action="store_true",
        help="If set, then compress the output directory into a Zip file.",
    )

    return parser.parse_args()


def clear_data(args):
    """Clear data path.

    If file_loc is default, then remove directly. Otherwise, prompt for
    acknowledgement.

    """
    if args.clear_data:
        if args.file_loc is None:
            shutil.rmtree(GenePlexus(log_level="CRITICAL").file_loc)
        else:
            file_loc = normexpand(args.file_loc)
            if input("Remove directory {file_loc}? [y/n]") == "y":
                shutil.rmtree(file_loc)
        exit()


def main():
    """Run the full GenePlexus pipeline."""
    args = parse_args()
    log_level = "CRITICAL" if args.quiet else args.log_level

    clear_data(args)  # data cleared if args.clear_data is true

    if "," in args.sp_res:
        args.sp_res = args.sp_res.split(",")
    if "," in args.gsc_res:
        args.gsc_res = args.gsc_res.split(",")

    # Create geneplexus object and auto download data files
    gp = GenePlexus(
        file_loc=args.file_loc,
        net_type=args.network,
        features=args.feature,
        sp_trn=args.sp_trn,
        sp_res=args.sp_res,
        gsc_trn=args.gsc_trn,
        gsc_res=args.gsc_res,
        auto_download=args.auto_download,
        log_level=log_level,
        log_to_file=True,
    )

    # Load input gene list
    gp.load_genes(read_gene_list(args.input_file, args.gene_list_delimiter))

    # load negative gene list if one provided
    if args.input_negatives != None:
        gp.load_negatives(read_gene_list(args.input_negatives, args.gene_list_delimiter))

    # run the pipeline
    if args.do_clustering:
        gp.cluster_input(
            clust_min_size=args.clust_min_size,
            clust_max_size=args.clust_max_size,
            clust_max_tries=args.clust_max_tries,
            clust_res=args.clust_res,
            clust_weighted=args.clust_unweighted,
        )
    gp.fit(
        logreg_kwargs=args.fit_logreg_kwargs,
        scale=args.fit_scale,
        min_num_pos=args.fit_min_num_pos,
        min_num_pos_cv=args.fit_min_num_pos_cv,
        num_folds=args.fit_num_folds,
        null_val=args.fit_null_val,
        random_state=args.fit_random_state,
        cross_validate=args.fit_skip_cross_validate,
    )
    gp.predict()
    if not args.skip_mdl_sim:
        gp.make_sim_dfs()
    else:
        logger.info("Skipping model similarity computation.")
    if not args.skip_sm_edgelist:
        gp.make_small_edgelist(
            num_nodes=args.small_edgelist_num_nodes,
        )
    else:
        logger.info("Skipping making small edgelist.")
    gp.save_class(
        args.output_dir,
        save_type=args.save_type,
        zip_output=args.zip_output,
        overwrite=args.overwrite,
    )
    gp.remove_log_file()


if __name__ == "__main__":
    main()
