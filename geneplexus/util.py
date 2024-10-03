"""Utilities including file and path handling."""
import functools
import json
import os
import os.path as osp
from threading import Thread
from typing import Any
from typing import Dict
from typing import Generator
from typing import List
from typing import Literal
from typing import Optional

import numpy as np

from . import config


def timeout(timeout: int, msg: str = ""):
    """Timeout decorator using thread join timeout.

    Args:
        timeout: Max function execution time in seconds.

    """

    def decorate(func, /):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            res = [TimeoutError(f"({timeout=}) {msg}")]

            def wraped_func():
                try:
                    res[0] = func(*args, **kwargs)
                except Exception as e:
                    res[0] = e

            t = Thread(target=wraped_func)
            t.daemon = True
            t.start()
            t.join(timeout)

            ret = res[0]
            if isinstance(ret, BaseException):
                print("")
                raise ret

            return ret

        return wrapper

    return decorate


def normexpand(path: str, create: bool = True) -> str:
    """Normalize then expand path and optionally create dir."""
    new_path = osp.abspath(osp.normpath(osp.expanduser(path)))
    if create:
        os.makedirs(new_path, exist_ok=True)
    return new_path


def format_choices(choices: List[str]) -> str:
    """Convert list of str to choices format."""
    return f"The choices are: {{{', '.join(choices)}}}"


def mapgene(gene: str, entrez_to_other: Dict[str, List[str]]) -> str:
    """Map entrez to other representations.

    Args:
        gene: Entrez gene ID.
        entrez_to_other: Mapping from Entrez to list of
            other gene representations of interest.

    Returns:
        str: Gene representation corresponding to the gene Entrez ID.

    Note:
        Mapping from a single Entrez to multiple representations is allowed and
        the representations will be separated by '/'.

    """
    try:
        syms = "/".join(entrez_to_other[gene])
    except KeyError:
        syms = "N/A"
    return syms


def get_all_filenames() -> Generator[str, None, None]:
    """Iterate over filenames."""
    with open(config.DATA_FILENAMES_PATH) as f:
        for line in f:
            yield line.strip()


def check_file(path: str):
    """Check existence of a file.

    Args:
        path: Path to the file.

    Raises:
        FileNotFoundError: if file not exist.

    """
    if not osp.isfile(path):
        raise FileNotFoundError(path)


def read_gene_list(
    path: str,
    sep: Optional[str] = "newline",
) -> List[str]:
    """Read gene list from flie.

    Args:
        path: Path to the input gene list file.
        sep: Seperator between genes (default: "newline").

    """
    if sep == "newline":
        sep = None
    elif sep == "tab":
        sep = "\t"
    return [gene.strip("'") for gene in open(path).read().split(sep)]


def _load_json_file(file_loc: str, file_name: str) -> Dict[str, Any]:
    """Load JSON into dictionary.

    Args:
        file_loc: Location of data files.
        file_name: Name of the file.

    """
    file_path = osp.join(file_loc, file_name)
    check_file(file_path)
    return json.load(open(file_path, "rb"))


def load_geneid_conversion(
    file_loc: str,
    species: config.SPECIES_TYPE,
    src_id_type: config.ID_SRC_TYPE,
    dst_id_type: config.ID_DST_TYPE,
    upper: bool = False,
) -> config.ID_CONVERSION_MAP_TYPE:
    """Load the gene ID conversion mapping.

    Args:
        file_loc: Directory containig the ID conversion file.
        species: The species of the ID conversion file.
        src_id_type: Souce gene ID type.
        dst_id_type: Destination gene ID type.
        upper: If set to True, then convert all keys to upper case.

    """
    if (src_id_type, dst_id_type) not in config.VALID_ID_CONVERSION:
        raise ValueError(f"Invalid ID conversion from {src_id_type} to {dst_id_type}")

    file_name = f"IDconversion__{species}__{src_id_type}-to-{dst_id_type}.json"
    conversion_map = _load_json_file(file_loc, file_name)

    if upper:
        conversion_map = {src.upper(): dst for src, dst in conversion_map.items()}

    return conversion_map


def load_gsc(
    file_loc: str,
    species: config.SPECIES_TYPE,
    gsc: config.GSC_TYPE,
    net_type: config.NET_TYPE,
) -> config.GSC_DATA_TYPE:
    """Load gene set collection dictionary.

    Args:
        file_loc: Location of data files.
        species: The species of the gene set collection.
        gsct: Target gene set collection.
        net_type: Network used.

    """
    file_name = f"GSC__{species}__{gsc}__{net_type}.json"
    return _load_json_file(file_loc, file_name)


def load_biomart(
    file_loc: str,
    sp_trn: config.SPECIES_TYPE,
    sp_tst: config.SPECIES_TYPE,
) -> config.BIOMART_DATA_TYPE:
    """Load biomart dictionary.

    Args:
        file_loc: Location of data files.
        sp_trn: The species used for training.
        sp_tst: The species in which the results are in.

    """
    file_name = f"BioMart__{sp_trn}__{sp_tst}.json"
    return _load_json_file(file_loc, file_name)


def load_pretrained_weights(
    file_loc: str,
    species: config.SPECIES_TYPE,
    gsc: config.GSC_TYPE,
    net_type: config.NET_TYPE,
    features: config.FEATURE_TYPE,
) -> config.PRETRAINED_DATA_TYPE:
    """Load pretrained model dictionary.

    Args:
        file_loc: Location of data files.
        species: The species of the gene set collection.
        gsc: The gene set collection.
        net_type: Network used.
        features: Type of features used.

    """
    file_name = f"PreTrainedWeights__{species}__{gsc}__{net_type}__{features}.json"
    return _load_json_file(file_loc, file_name)


def _load_np_file(
    file_loc: str,
    file_name: str,
    load_method: Literal["npy", "txt"],
) -> np.ndarray:
    """Check np file existence and load.

    Args:
        file_loc: Location of data files.
        file_name: Name of the file.
        load_method: How to load the file ('npy' or 'txt').

    """
    file_path = osp.join(file_loc, file_name)
    check_file(file_path)

    if load_method == "npy":
        return np.load(file_path)
    elif load_method == "txt":
        return np.loadtxt(file_path, dtype=str)
    else:
        raise ValueError(f"Unknwon load method: {load_method!r}")


def load_node_order(
    file_loc: str,
    species: config.SPECIES_TYPE,
    net_type: config.NET_TYPE,
) -> np.ndarray:
    """Load network genes.

    Args:
        file_loc: Location of data files.
        species: The species of files.
        net_type: Network used.

    """
    file_name = f"NodeOrder__{species}__{net_type}.txt"
    return _load_np_file(file_loc, file_name, load_method="txt")


def load_gene_features(
    file_loc: str,
    species: config.SPECIES_TYPE,
    features: config.FEATURE_TYPE,
    net_type: config.NET_TYPE,
) -> np.ndarray:
    """Load gene features.

    Args:
        file_loc: Location of data files.
        species: The species of files.
        features: Type of features used.
        net_type: Network used.

    """
    file_name = f"Data__{species}__{features}__{net_type}.npy"
    return _load_np_file(file_loc, file_name, load_method="npy")


def param_warning(
    param_type: str,
    user_param: str,
    param_options: List[str],
) -> str:
    """Make custom param warning.

    Args:
        param_type: The parameter type.
        user_param: The user supplied parameter.
        param_options: List of options package natively supports.

    """
    return (
        f"\n{user_param} appears to be a custom {param_type} type. "
        f"PyGenePlexus natively supported {param_type} options are "
        f"{param_options}.Please make sure all custom files are named "
        "and formatted correctly. See docs for more information"
        "on custom files. "
    )


def get_all_net_types(file_loc: Optional[str], species: str) -> List[str]:
    """Return list of networks found in the data directory.

    Note:
        Only the node ordering files are checked (starts with ``NodeOrder``).

    """
    all_net_types = set()
    if file_loc:
        all_net_types.update(
            [i.split("__")[2].split(".txt")[0] for i in os.listdir(file_loc) if i.startswith(f"NodeOrder__{species}")],
        )
    return sorted(all_net_types)
