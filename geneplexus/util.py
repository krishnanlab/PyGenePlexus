import os.path as osp
import pickle
from typing import Dict
from typing import List

from . import config


def check_file(path: str):
    """Check existence of a file.

    Args:
        path (str): Path to the file.

    Raise:
        FileNotFoundError: if file not exist.

    """
    if not osp.isfile(path):
        raise FileNotFoundError(path)


def read_gene_list(
    path: str,
    sep: str = ", ",
) -> List[str]:
    """Read gene list from flie.

    Args:
        path (str): Path to the input gene list file.
        sep (str): Seperator between genes.

    """
    return [gene.strip("'") for gene in open(path, "r").read().split(sep)]


def get_geneid_conversion(
    file_loc: str,
    src_id_type: config.ID_SRC_TYPE,
    dst_id_type: config.ID_DST_TYPE,
    upper: bool = False,
) -> Dict[str, List[str]]:
    """Obtain the gene ID conversion mapping.

    Args:
        file_loc (str): Directory containig the ID conversion file.
        src_id_type (ID_SRC_TYPE): Souce gene ID type.
        dst_id_type (ID_DST_TYPE): Destination gene ID type.
        upper (bool): If set to True, then convert all keys to upper case.

    """
    if (src_id_type, dst_id_type) not in config.VALID_ID_CONVERSION:
        raise ValueError(f"Invalid ID conversion from {src_id_type} to {dst_id_type}")

    file_name = f"IDconversion_Homo-sapiens_{src_id_type}-to-{dst_id_type}.pickle"
    file_path = osp.join(file_loc, file_name)
    check_file(file_path)

    with open(file_path, "rb") as handle:
        conversion_map = pickle.load(handle)

    if upper:
        conversion_map = {src.upper(): dst for src, dst in conversion_map.items()}

    return conversion_map
