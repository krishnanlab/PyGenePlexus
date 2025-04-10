"""Command line interface for the GenePlexus pipeline."""
import argparse
import atexit
import json
import os
import os.path as osp
import pathlib
import shutil
import tempfile

import numpy as np
import pandas as pd

from . import config
from ._config import logger
from ._config.logger_util import attach_file_handler
from .geneplexus import GenePlexus
from .util import format_choices
from .util import normexpand
from .util import read_gene_list

os.environ["COLUMNS"] = "100"  # for CLI help page wrap line

TMP_LOG_FP, TMP_LOG_PATH = tempfile.mkstemp(suffix="_run.log")
FILE_HANDLER = attach_file_handler(logger, log_path=TMP_LOG_PATH)


def parse_args() -> argparse.Namespace:
    """Parse arguments from command line."""
    parser = argparse.ArgumentParser(
        description="Run the GenePlexus pipline on a input gene list.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=argparse.SUPPRESS,
    )

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
        "tabs. Other generic separator are also supported, e.g. ', '.",
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
        "-s",
        "--small_edgelist_num_nodes",
        default=50,
        metavar="",
        type=int,
        help="Number of nodes in the small edgelist.",
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
        "-ad",
        "--auto_download_off",
        action="store_true",
        help="Turns off autodownloader which is on by default.",
    )

    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress log messages (same as setting log_level to CRITICAL).",
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
        help="Overwrite existing result directory if set.",
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

    parser.add_argument(
        "--do_clustering",
        action="store_true",
        help="Do clustering step.",
    )

    parser.add_argument(
        "-st",
        "--save_type",
        default="all",
        metavar="",
        type=str,
        help="Which files to save.",
    )

    return parser.parse_args()


def run_pipeline(
    gp: GenePlexus,
    do_clustering: bool,
    num_nodes: int,
    skip_mdl_sim: bool,
    skip_sm_edgelist: bool,
    output_dir: str,
    save_type: str,
    zip_output: bool,
    overwrite: bool,
):
    """Run the full GenePlexus pipeline.

    Args:
        do_clustering: wehter to run the clustering step or not
        num_nodes: Number of top predicted genes to include in the induced
            subgraph.
        skip_mdl_sim: Whether or not to skip the computation of model
            similarities with GO and Mondo. This option is not yet available
            for custom networks.
        skip_sm_edgelist: skips making small edgelist of top predictions

    """
    if do_clustering:
        gp.cluster_input()
    gp.fit()
    gp.predict()
    if not skip_mdl_sim:
        gp.make_sim_dfs()
    else:
        logger.info("Skipping model similarity computation.")
    if not skip_sm_edgelist:
        gp.make_small_edgelist(num_nodes=num_nodes)
    else:
        logger.info("Skipping making small edgelist.")
    gp.save_class(
        output_dir,
        save_type=save_type,
        zip_output=zip_output,
        overwrite=overwrite,
    )


# def df_to_tsv(df: pd.DataFrame, root: str, name: str):
#     """Save a dataframe as a tsv file.
#
#     Args:
#         df: DataFrame to be saved.
#         root: Output directory,
#         name: Name of the file to be saved.
#
#     """
#     df.to_csv(osp.join(root, name), sep="\t", index=False)


# def save_results(gp, outdir, zip_output, overwrite, skip_mdl_sim, skip_sm_edgelist):
#     """Save all results generated by the GenePlexus pipeline.
#
#     Args:
#         outdir: Output directory.
#         zip_output: Whether or not to zip the output directory into a zip file.
#         overwrite: Whether or not to overwrite existing results.
#         skip_mdl_sim: Whether or not to skip the computation of model
#             similarities with GO, Monarch and/or Mondo.
#         skip_sm_edgelist: Whether or not to skip making the small edgelist.
#
#     """
#     zip_outpath = _suffix_fn(f"{outdir}.zip", overwrite=overwrite)
#     outdir = _suffix_dir(outdir, overwrite=overwrite, mktmp=zip_output)
#     df_to_tsv(gp.df_convert_out, outdir, "df_convert_out.tsv")
#     np.savetxt(osp.join(outdir, "pos_genes_in_net.txt"), gp.pos_genes_in_net, fmt="%s")
#     np.savetxt(osp.join(outdir, "negative_genes.txt"), gp.negative_genes, fmt="%s")
#     np.savetxt(osp.join(outdir, "net_genes.txt"), gp.net_genes, fmt="%s")
#     with open(osp.join(outdir, "neutral_gene_info.json"), "w") as f:
#         json.dump(gp.neutral_gene_info, f)
#     np.savetxt(osp.join(outdir, "avgps.txt"), gp.avgps, fmt="%.18f")
#     np.savetxt(osp.join(outdir, "mdl_weights.txt"), gp.mdl_weights, fmt="%.18f")
#     df_to_tsv(gp.df_probs, outdir, "df_probs.tsv")
#     if not skip_mdl_sim:
#         df_to_tsv(gp.df_sim, outdir, "df_sim.tsv")
#     if not skip_mdl_sim:
#         df_to_tsv(gp.df_edge, outdir, "df_edge.tsv")
#         df_to_tsv(gp.df_edge_sym, outdir, "df_edge_sym.tsv")
#         np.savetxt(osp.join(outdir, "isolated_genes.txt"), gp.isolated_genes, fmt="%s")
#         np.savetxt(osp.join(outdir, "isolated_genes_sym.txt"), gp.isolated_genes_sym, fmt="%s")
#     df_to_tsv(gp.df_convert_out_subset, outdir, "df_convert_out_subset.tsv")
#
#     # Dump config, close file handler and move run log to result directory
#     gp.dump_config(outdir)
#     logger.removeHandler(FILE_HANDLER)
#     FILE_HANDLER.flush()
#     FILE_HANDLER.close()
#     os.close(TMP_LOG_FP)  # https://stackoverflow.com/a/60357401
#     shutil.move(TMP_LOG_PATH, osp.join(outdir, "run.log"))
#
#     # Optionally zip the result directory
#     if zip_output:
#         outpath = pathlib.Path(outdir)
#         logger.info("Zipping output files")
#         shutil.make_archive(zip_outpath[:-4], "zip", outpath.parent, outpath.name)
#         shutil.rmtree(outdir)
#         logger.info(f"Removing temporary directory {outdir}")
#         logger.info(f"Done! Results saved to {zip_outpath}")
#     else:
#         logger.info(f"Done! Results saved to {outdir}")


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


# def _suffix_dir(path, idx=0, overwrite=False, mktmp=False):
#     """Add int suffix to dir name if nonempty dir existed."""
#     if mktmp:
#         new_path = tempfile.mkdtemp(dir=pathlib.Path(path).absolute().parent)
#         logger.info(
#             f"Results to be zipped, create temporary directory for holding unzipped results {new_path}",
#         )
#         return new_path
#     new_path = normexpand(f"{path}_{idx}" if idx > 0 else path)
#     if os.listdir(new_path):
#         if overwrite:
#             logger.warning(f"Output directory exits {path}, overwriting.")
#             shutil.rmtree(new_path)
#             os.makedirs(new_path)
#         else:
#             new_path = _suffix_dir(path, idx=idx + 1)
#     elif path != new_path:
#         logger.warning(f"Output directory exists {path}, redirecting to {new_path}")
#     return new_path


# def _suffix_fn(path, idx=0, overwrite=False):
#     """Add int suffix to file name if file existed."""
#     new_path = f"_{idx}".join(osp.splitext(path)) if idx > 0 else path
#     if osp.isfile(new_path):
#         if overwrite:
#             logger.warning(f"Output zip file exits {path}, overwriting.")
#             os.remove(new_path)
#         else:
#             new_path = _suffix_fn(path, idx=idx + 1)
#     elif path != new_path:
#         logger.warning(f"Output zip file exists {path}, redirecting to {new_path}")
#     return new_path


@atexit.register
def interrupted():
    """Check if program is interrupted and print temporary log file path."""
    if osp.isfile(TMP_LOG_PATH):
        logger.critical(f"Program interrupted, temporary run log saved at {TMP_LOG_PATH}")


def main():
    """Command line interface."""
    args = parse_args()
    log_level = "CRITICAL" if args.quiet else args.log_level
    if args.auto_download_off:
        auto_download = False
    else:
        auto_download = True

    clear_data(args)  # data cleared if args.clear_data is true

    if "," in args.sp_res:
        args.sp_res = args.sp_res.split(",")
    if "," in args.gsc_res:
        args.gsc_res = args.gsc_res.split(",")

    # Create geneplexus object and auto download data files
    gp = GenePlexus(
        file_loc=args.data_dir,
        net_type=args.network,
        features=args.feature,
        sp_trn=args.sp_trn,
        sp_res=args.sp_res,
        gsc_trn=args.gsc_trn,
        gsc_res=args.gsc_res,
        auto_download=auto_download,
        log_level=log_level,
    )

    # Load input gene list
    gp.load_genes(read_gene_list(args.input_file, args.gene_list_delimiter))

    # Run pipeline and save results
    run_pipeline(
        gp,
        args.do_clustering,
        args.small_edgelist_num_nodes,
        args.skip_mdl_sim,
        args.skip_sm_edgelist,
        args.output_dir,
        args.save_type,
        args.zip_output,
        args.overwrite,
    )
    # save_results(
    #     gp,
    #     normexpand(args.output_dir),
    #     args.zip_output,
    #     args.overwrite,
    #     args.skip_mdl_sim,
    #     args.skip_sm_edgelist,
    # )


if __name__ == "__main__":
    main()
