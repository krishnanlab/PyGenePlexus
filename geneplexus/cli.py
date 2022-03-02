"""Command line interface for the GenePlexus pipeline."""
import argparse
import os
import os.path as osp
import pathlib
from typing import Tuple

import pandas as pd

from . import config
from .download import download_select_data
from .geneplexus import GenePlexus
from .util import read_gene_list


HOMEDIR = pathlib.Path(__file__).absolute().parent.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the GenePlexus pipline on a input gene list.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_file",
        required=True,
        help="Input gene list (.txt) file (one gene per line).",
    )

    parser.add_argument(
        "--network",
        default="BioGRID",
        choices=config.ALL_NETWORKS,
        help="Network to use for generating features.",
    )

    parser.add_argument(
        "--feature",
        default="Embedding",
        choices=config.ALL_FEATURES,
        help="Types of feature to use.",
    )

    parser.add_argument(
        "--GSC",
        default="GO",
        choices=config.ALL_GSCS,
        help="Geneset collection used to generate negatives and the model similarities.",
    )

    parser.add_argument(
        "--small_edgelist_num_nodes",
        default=50,
        type=int,
        help="Number of nodes in the small edgelist.",
    )

    parser.add_argument(
        "--data_dir",
        default="data/",
        help="Directory in which the data are stored.",
    )

    parser.add_argument(
        "--output_dir",
        default="result/",
        help="Output directory with respect to the repo root directory.",
    )

    return parser.parse_args()


def preprocess(args: argparse.Namespace) -> Tuple[str, str]:
    datadir = osp.join(HOMEDIR, args.data_dir)
    outdir = osp.join(HOMEDIR, args.output_dir)
    os.makedirs(datadir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)
    download_select_data(datadir, "All", args.network, args.feature, ["GO", "DisGeNet"])
    return datadir, outdir


def run_pipeline(gp: GenePlexus, num_nodes: int):
    gp.convert_to_Entrez()
    gp.get_pos_and_neg_genes()
    gp.fit_and_predict()
    gp.make_sim_dfs()
    gp.make_small_edgelist(num_nodes=num_nodes)
    gp.alter_validation_df()


def df_to_tsv(df: pd.DataFrame, root: str, name: str):
    df.to_csv(osp.join(root, name), sep="\t", index=False)


def save_results(gp, outdir):
    df_to_tsv(gp.df_convert_out, outdir, "df_convert_out.tsv")
    df_to_tsv(gp.df_probs, outdir, "df_probs.tsv")
    df_to_tsv(gp.df_sim_GO, outdir, "df_sim_GO.tsv")
    df_to_tsv(gp.df_sim_Dis, outdir, "df_sim_Dis.tsv")
    df_to_tsv(gp.df_edge, outdir, "df_edge.tsv")
    df_to_tsv(gp.df_edge_sym, outdir, "df_edge_sym.tsv")
    df_to_tsv(gp.df_convert_out_subset, outdir, "df_convert_out_subset.tsv")


def main():
    args = parse_args()
    datadir, outdir = preprocess(args)

    gp = GenePlexus(datadir, args.network, args.feature, args.GSC)
    gp.load_genes(read_gene_list(args.input_file))

    run_pipeline(gp, args.input_file)
    save_results(gp, outdir)


if __name__ == "__main__":
    main()
