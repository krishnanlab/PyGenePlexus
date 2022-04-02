"""Command line interface for the GenePlexus pipeline."""
import argparse
import logging
import os.path as osp
import pathlib
import shutil

import pandas as pd

from . import config
from ._config import logger
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
        default=None,
        metavar="",
        help="Directory in which the data are stored, if set to None, then use "
        "the default data directory ~/.data/geneplexus",
    )

    parser.add_argument(
        "-od",
        "--output_dir",
        default="result/",
        metavar="",
        help="Output directory with respect to the repo root directory.",
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

    parser.add_argument(
        "-z",
        "--zip-output",
        action="store_true",
        help="If set, then compress the output directory into a Zip file.",
    )

    parser.add_argument(
        "--clear-data",
        action="store_true",
        help="Clear data directory and exit.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing results under the specified result output "
        "directory when set. Otherwise, prompt for overwriting acknowledgement "
        "when the result file already exist.",
    )

    return parser.parse_args()


def run_pipeline(gp: GenePlexus, num_nodes: int):
    """Run the full GenePlexus pipeline."""
    gp.convert_to_Entrez()
    gp.get_pos_and_neg_genes()
    gp.fit_and_predict()
    gp.make_sim_dfs()
    gp.make_small_edgelist(num_nodes=num_nodes)
    gp.alter_validation_df()


def _ack_overwrite(outpath: str, overwrite: bool) -> bool:
    if osp.isfile(outpath):
        if overwrite:
            logger.warning(f"Overwritting file {outpath}")
        else:
            while True:
                ans = input(
                    f"File exists {outpath}, overwrite anyway? (specify the "
                    "--overwrite cli option to overwrite without prompt) [y/n]",
                )
                if ans == "y":
                    logger.info(f"Overwritting file {outpath}")
                    break
                elif ans == "n":
                    logger.info(f"Refected file overwriting {outpath}.")
                    return False
                else:
                    print("Please answer 'y' or 'n'.")
    return True


def df_to_tsv(df: pd.DataFrame, root: str, name: str, overwrite: bool):
    """Save a dataframe as a tsv file.

    Args:
        df: DataFrame to be saved.
        root: Output directory,
        name: Name of the file to be saved.
        overwrite: If set to True, overwrite file even if it already exists.
            Otherwise, prompt for overwriting acknowledgement.

    """
    outpath = osp.join(root, name)
    if _ack_overwrite(outpath, overwrite):
        df.to_csv(outpath, sep="\t", index=False)


def save_results(gp, outdir, overwrite, zip_output):
    """Save all results generated by the GenePlexus pipeline."""
    df_to_tsv(gp.df_convert_out, outdir, "df_convert_out.tsv", overwrite)
    df_to_tsv(gp.df_probs, outdir, "df_probs.tsv", overwrite)
    df_to_tsv(gp.df_sim_GO, outdir, "df_sim_GO.tsv", overwrite)
    df_to_tsv(gp.df_sim_Dis, outdir, "df_sim_Dis.tsv", overwrite)
    df_to_tsv(gp.df_edge, outdir, "df_edge.tsv", overwrite)
    df_to_tsv(gp.df_edge_sym, outdir, "df_edge_sym.tsv", overwrite)
    df_to_tsv(gp.df_convert_out_subset, outdir, "df_convert_out_subset.tsv", overwrite)

    # Optionally zip the result directory
    if zip_output:
        zip_outpath = f"{outdir}.zip"
        if _ack_overwrite(zip_outpath, overwrite):
            logger.info("Zipping output files")
            outpath = pathlib.Path(outdir)
            shutil.make_archive(outdir, "zip", outpath.parent, outpath.name)
            shutil.rmtree(outdir)
            logger.info(f"Done! Results saved to {zip_outpath}")
    else:
        logger.info(f"Done! Results saved to {outdir}")


def clear_data(args):
    """Clear data path.

    If data_dir is default, then remove directly. Otherwise, prompt for
    acknowledgement.

    """
    if args.clear_data:
        if args.data_dir is None:
            shutil.rmtree(GenePlexus(log_level="CRITICAL").file_loc)
        else:
            data_dir = normexpand(args.data_dir)
            if input("Remove directory {data_dir}? [y/n]") == "y":
                shutil.rmtree(data_dir)
        exit()


def main():
    """Command line interface."""
    args = parse_args()
    log_level = logging.getLevelName("CRITICAL" if args.quiet else args.log_level)

    clear_data(args)

    # Create geneplexus object and auto download data files
    gp = GenePlexus(
        args.data_dir,
        args.network,
        args.feature,
        args.gsc,
        auto_download=True,
        log_level=log_level,
    )

    # Load input gene list
    gp.load_genes(read_gene_list(args.input_file, args.gene_list_delimiter))

    # Run pipeline and save results
    run_pipeline(gp, args.input_file)
    save_results(gp, normexpand(args.output_dir), args.overwrite, args.zip_output)


if __name__ == "__main__":
    main()
