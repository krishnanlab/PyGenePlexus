"""Command line interface for the GenePlexus pipeline."""
import argparse
import logging
import os
import os.path as osp
import pathlib
import tarfile
from typing import Tuple

import pandas as pd

from . import config
from ._config import logger
from .download import download_select_data
from .geneplexus import GenePlexus
from .util import format_choices
from .util import normexpand
from .util import read_gene_list


def parse_args() -> argparse.Namespace:
    """Parse arguments from command line."""
    parser = argparse.ArgumentParser(
        description="Run the GenePlexus pipline on a input gene list.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input_file",
        metavar="",
        required=True,
        help="Input gene list (.txt) file (one gene per line).",
    )

    parser.add_argument(
        "-d",
        "--gene_list_delimiter",
        default="newline",
        metavar="",
        help="Delimiter used in the gene list. Use 'newline' if the genes are "
        "separated by new line, and use 'tab' if the genes are seperate by "
        "tabs. Other generic separator are also supported, e.g. ', '.",
    )

    parser.add_argument(
        "-n",
        "--network",
        default="BioGRID",
        metavar="",
        help="Network to use. {format_choices(config.ALL_NETWORKS)}",
    )

    parser.add_argument(
        "-f",
        "--feature",
        default="Embedding",
        metavar="",
        choices=config.ALL_FEATURES,
        help=f"Types of feature to use. {format_choices(config.ALL_FEATURES)}",
    )

    parser.add_argument(
        "-g",
        "--gsc",
        default="GO",
        metavar="",
        help="Geneset collection used to generate negatives and the model"
        f"similarities. {format_choices(config.ALL_GSCS)}",
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
        "-dd",
        "--data_dir",
        default="data/",
        metavar="",
        help="Directory in which the data are stored.",
    )

    parser.add_argument(
        "-od",
        "--output_dir",
        default="result/",
        metavar="",
        help="Output directory with respect to the repo root directory.",
    )

    parser.add_argument(
        "-z",
        "--zip_output",
        action="store_true",
        help="If set, then compress the output directory into a Tar Gz file.",
    )

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
        help="Suppress log messages (same as setting lov_level to CRITICAL).",
    )

    return parser.parse_args()


def preprocess(args: argparse.Namespace) -> Tuple[str, str]:
    """Set up data and result directories and download data if necessary."""
    datadir = normexpand(args.data_dir)
    outdir = normexpand(args.output_dir)
    os.makedirs(datadir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)
    download_select_data(datadir, "All", args.network, args.feature, ["GO", "DisGeNet"])
    return datadir, outdir


def run_pipeline(gp: GenePlexus, num_nodes: int):
    """Run the full GenePlexus pipeline."""
    gp.convert_to_Entrez()
    gp.get_pos_and_neg_genes()
    gp.fit_and_predict()
    gp.make_sim_dfs()
    gp.make_small_edgelist(num_nodes=num_nodes)
    gp.alter_validation_df()


def df_to_tsv(df: pd.DataFrame, root: str, name: str):
    """Save a dataframe as a tsv file."""
    df.to_csv(osp.join(root, name), sep="\t", index=False)


def save_results(gp, outdir, zip_output):
    """Save all results generated by the GenePlexus pipeline."""
    df_to_tsv(gp.df_convert_out, outdir, "df_convert_out.tsv")
    df_to_tsv(gp.df_probs, outdir, "df_probs.tsv")
    df_to_tsv(gp.df_sim_GO, outdir, "df_sim_GO.tsv")
    df_to_tsv(gp.df_sim_Dis, outdir, "df_sim_Dis.tsv")
    df_to_tsv(gp.df_edge, outdir, "df_edge.tsv")
    df_to_tsv(gp.df_edge_sym, outdir, "df_edge_sym.tsv")
    df_to_tsv(gp.df_convert_out_subset, outdir, "df_convert_out_subset.tsv")

    # Optionally zip the result directory
    if zip_output:
        file_name = pathlib.Path(outdir).name
        with tarfile.open(f"{file_name}.tar.gz", "w:gz") as tar:
            tar.add(pathlib.Path(outdir), arcname=file_name)


def main():
    """Command line interface."""
    args = parse_args()
    logger.setLevel(logging.getLevelName("CRITICAL" if args.quiet else args.log_level))
    datadir, outdir = preprocess(args)

    gp = GenePlexus(datadir, args.network, args.feature, args.gsc)
    gp.load_genes(read_gene_list(args.input_file, args.gene_list_delimiter))

    run_pipeline(gp, args.input_file)
    save_results(gp, outdir, args.zip_output)


if __name__ == "__main__":
    main()
