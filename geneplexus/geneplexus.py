"""GenePlexus API."""
import os
import os.path as osp
import tempfile
import warnings
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

import numpy as np
import pystow
import yaml

from . import _geneplexus
from . import util
from ._config import config
from ._config import logger
from ._config.logger_util import attach_file_handler
from ._config.logger_util import set_stream_level
from .download import download_select_data
from .exception import NoPositivesError


class ModelInfo:
    def __init__(self):
        """Class to hold the trainig objects"""
        pass


class ModelResults:
    def __init__(self):
        """Class to hold the result objects"""
        pass


class GenePlexus:
    """The GenePlexus API class."""

    def __init__(
        self,
        file_loc: Optional[str] = config.DEFAULT_PARAMETERS["file_loc"],
        net_type: config.NET_TYPE = config.DEFAULT_PARAMETERS["net_type"],
        features: config.FEATURE_TYPE = config.DEFAULT_PARAMETERS["features"],
        sp_trn: config.SPECIES_TYPE = config.DEFAULT_PARAMETERS["sp_trn"],
        sp_res: config.SPECIES_SELECTION_TYPE = config.DEFAULT_PARAMETERS["sp_res"],
        gsc_trn: config.GSC_TYPE = config.DEFAULT_PARAMETERS["gsc_trn"],
        gsc_res: config.GSC_SELECTION_TYPE = config.DEFAULT_PARAMETERS["gsc_res"],
        input_genes: Optional[List[str]] = config.DEFAULT_PARAMETERS["input_genes"],
        input_negatives: Optional[List[str]] = config.DEFAULT_PARAMETERS["input_negatives"],
        auto_download: bool = config.DEFAULT_PARAMETERS["auto_download"],
        log_level: config.LOG_LEVEL_TYPE = config.DEFAULT_PARAMETERS["log_level"],
        log_to_file: bool = config.DEFAULT_PARAMETERS["log_to_file"],
    ):
        """Initialize the GenePlexus object.

        Args:
            file_loc: Location of data files, if not specified, set to default
                data path ``~/.data/geneplexus``
            net_type: Type of network to use
            features: Type of features of the network to use
            sp_trn: The species of the training data
            sp_res: The species the results are in, can be a list
            gsc_trn: Gene set collection used during training
            gsc_res: Gene set(s) collection used when generating
                results, can be a list. If list needs to be the same
                length as number of results species to do
            input_genes: Input gene list, can be mixed type. Can also be set
                later if not specified at init time by simply calling
                :meth:`load_genes`.
            input_negatives: Input list of negative genes, can be mixed type.
                Can also be set later if not specified at init time by simply calling
                :meth:`load_negatives`.
            auto_download: Automatically download necessary files if set.
            log_level: Logging level.
            log_to_file: If True logger will be saved when `save_class` is run. This
                will also create a tmp file that can be explicitly deleted with
                `remove_log_file`.

        The following clsss attributes are set when ``__init__`` is run

        :attr:`GenePlexus._is_custom` bool
            If the species, network or feature type was supplied by the user.
        :attr:`GenePlexus._file_loc` str
            File path set for the data.
        :attr:`GenePlexus._features` str
            Type of network features used.
        :attr:`GenePlexus._sp_trn` str
            Species used in training.
        :attr:`GenePlexus._sp_res` (str, List[str])
            Species used in the results.
        :attr:`GenePlexus._gsc_trn` str
            Gene set collection used in training.
        :attr:`GenePlexus._gsc_res` (str, List[str])
            TGene set collection(s) used in results.
        :attr:`GenePlexus._net_type` str
            Type of network used.
        :attr:`GenePlexus.log_level` str
            The verbosity of the logger.
        :attr:`GenePlexus.log_to_file` bool
            Whether or not the log file was saved as a file
        :attr:`GenePlexus.auto_download` bool
            If data was attmepted to be auto downloaded.
        :attr:`GenePlexus.gsc_trn_original` str
            If internal data checks are run, this can different than _gsc_trn.
        :attr:`GenePlexus.gsc_res_original` (str, List[str])
            If internal data checks are run, this can different that _gsc_res.
        :attr:`GenePlexus.sp_gsc_pairs` List[str]
            The combination of all sp and gsc used, hyphen separated.
        :attr:`GenePlexus.model_info[ModelName]` Class
            model_info is a dictionary where each key is a different model and holds the ModelInfo class.
        :attr:`GenePlexus.model_info[ModelName].results[ResultName]` Class
            results is a dictionary where each key is a different result and holds ModelResults class.


        """
        set_stream_level(logger, log_level)
        if log_to_file:
            # Create a temporary file that is deleted on close or script exit
            TMP_LOG_FP, TMP_LOG_PATH = tempfile.mkstemp(suffix="_geneplexus.log")
            logger.info(f"The TMP_LOG_PATH is {TMP_LOG_PATH}")
            FILE_HANDLER = attach_file_handler(logger, log_path=TMP_LOG_PATH, log_level=log_level)
            self.file_handler = FILE_HANDLER
            self.log_tmp_path = TMP_LOG_PATH
        self._is_custom: bool = False
        self.file_loc = file_loc  # type: ignore
        self.features = features
        self.sp_trn = sp_trn
        self.sp_res = sp_res
        self.gsc_trn = gsc_trn
        self.gsc_res = gsc_res
        self.net_type = net_type
        self.log_level = log_level
        self.log_to_file = log_to_file
        self.auto_download = auto_download
        self.input_genes: List[str] = input_genes
        self.input_negatives: List[str] = input_negatives

        if self.auto_download and self._is_custom:
            warnings.warn(
                "\nSkipping auto download for custom files. Unset auto_download option to suppress this message.",
                UserWarning,
                stacklevel=2,
            )
        elif self.auto_download:
            download_select_data(
                self.file_loc,
                [self.sp_trn] + self.sp_res,
                log_level=log_level,
            )

        if self._is_custom:
            warnings.warn(
                f"is_custom is set to True either manually "
                "or by autodection of species or GSC not "
                "contained in the pre-processed data. All "
                "compatability checks are being turned off.",
            )
            self.gsc_trn_original = self.gsc_trn
            self.gsc_res_original = self.gsc_res
        else:
            # check option compatability for preprocessed data
            sp_res_subset, gsc_res_subset = util.data_checks(
                self.sp_trn,
                self.net_type,
                self.gsc_trn,
                self.sp_res,
                self.gsc_res,
            )
            self.sp_res = sp_res_subset
            self.gsc_res = gsc_res_subset

            # for combined, display contexts and change some GSC names
            gsc_trn_updated, gsc_res_updated = util.combined_info(
                self.sp_trn,
                self.gsc_trn,
                self.sp_res,
                self.gsc_res,
            )
            self.gsc_trn_original = self.gsc_trn
            self.gsc_trn = gsc_trn_updated
            self.gsc_res_original = self.gsc_res
            self.gsc_res = gsc_res_updated

        # remove duplicate sp-gsc combos if any for results
        sp_res_nodup, gsc_res_nodup, gsc_res_original_nodup = util.remove_duplicates(
            self.sp_res,
            self.gsc_res,
            self.gsc_res_original,
        )
        self.sp_res = sp_res_nodup
        self.gsc_res = gsc_res_nodup
        self.gsc_res_original = gsc_res_original_nodup

        # create objects that will always be used
        self.sp_gsc_pairs = ["-".join(str(item) for item in pair) for pair in zip(self.sp_res, self.gsc_res_original)]
        self.model_info = {"All-Genes": ModelInfo()}
        self.model_info["All-Genes"].results = {}
        for apair in self.sp_gsc_pairs:
            self.model_info["All-Genes"].results[apair] = ModelResults()

        # set a clus_min_size to make sure it matching min_num_pos later
        self.clust_min_size = None

        if input_genes is not None:
            self.load_genes(input_genes)

        if input_negatives is not None:
            self.load_negatives(input_negatives)

    @property
    def file_loc(self) -> str:
        """File location.
        Use default data location ~/.data/geneplexus if not set.
        """
        return self._file_loc

    @file_loc.setter
    def file_loc(self, file_loc: Optional[str]):
        if file_loc is None:
            self._file_loc = str(pystow.join("geneplexus"))
        else:
            self._file_loc = util.normexpand(file_loc)
        logger.info(f"Data direcory set to {self._file_loc}")

    @property
    def net_type(self) -> config.NET_TYPE:
        """Network to use."""
        return self._net_type

    @net_type.setter
    def net_type(self, net_type: config.NET_TYPE):
        if net_type not in config.ALL_NETWORKS:
            warnings.warn(
                util.param_warning("network", net_type, config.ALL_NETWORKS),
                UserWarning,
                stacklevel=2,
            )
            self._is_custom = True
            logger.info(f"Using custom network {net_type!r}")

        self._net_type = net_type

    @property
    def features(self) -> config.FEATURE_TYPE:
        """Features to use."""
        return self._features

    @features.setter
    def features(self, features: config.FEATURE_TYPE):
        if features not in config.ALL_FEATURES:
            warnings.warn(
                util.param_warning("feature", features, config.ALL_FEATURES),
                UserWarning,
                stacklevel=2,
            )
            self._is_custom = True
            logger.info(f"Using custom feature {features!r}")
        self._features = features

    @property
    def sp_trn(self) -> config.SPECIES_TYPE:
        """Training species."""
        return self._sp_trn

    @sp_trn.setter
    def sp_trn(self, sp_trn: config.SPECIES_TYPE):
        if sp_trn not in config.ALL_SPECIES:
            warnings.warn(
                util.param_warning("species", sp_trn, config.ALL_SPECIES),
                UserWarning,
                stacklevel=2,
            )
            self._is_custom = True
            logger.info(f"Using custom species {sp_trn!r}")
        self._sp_trn = sp_trn

    @property
    def sp_res(self) -> config.SPECIES_SELECTION_TYPE:
        """Results_species."""
        return self._sp_res

    @sp_res.setter
    def sp_res(self, sp_res: config.SPECIES_SELECTION_TYPE):
        if isinstance(sp_res, str):
            if sp_res == "All":
                sp_res = config.ALL_SPECIES
            else:
                sp_res = [sp_res]
        elif not isinstance(sp_res, list):
            raise TypeError(f"Expected str type or list of str type, got {type(sp_res)}")
        for i in sp_res:
            if i not in config.ALL_SPECIES:
                warnings.warn(
                    util.param_warning("species", i, config.ALL_SPECIES),
                    UserWarning,
                    stacklevel=2,
                )
                self._is_custom = True
                logger.info(f"There is a custom species in {sp_res!r}")
        self._sp_res = sp_res

    @property
    def gsc_trn(self) -> config.GSC_TYPE:
        """Geneset collection used in training."""
        return self._gsc_trn

    @gsc_trn.setter
    def gsc_trn(self, gsc_trn: config.GSC_TYPE):
        if gsc_trn not in config.ALL_GSCS:
            warnings.warn(
                util.param_warning("GSC", gsc_trn, config.ALL_GSCS),
                UserWarning,
                stacklevel=2,
            )
            self._is_custom = True
            logger.info(f"Using custom GSC {gsc_trn!r}")
        self._gsc_trn = gsc_trn

    @property
    def gsc_res(self) -> config.GSC_SELECTION_TYPE:
        """Geneset collection used when generating results."""
        return self._gsc_res

    @gsc_res.setter
    def gsc_res(self, gsc_res: config.GSC_SELECTION_TYPE):
        if isinstance(gsc_res, str):
            gsc_res = [gsc_res] * len(self.sp_res)
        elif not isinstance(gsc_res, list):
            raise TypeError(f"Expected str type or list of str type, got {type(gsc_res)}")
        if len(self.sp_res) != len(gsc_res):
            raise ValueError(f"Length of sp_res list not the same as gsc_res list")
        for i in gsc_res:
            if i not in config.ALL_GSCS:
                warnings.warn(
                    util.param_warning("GSC", i, config.ALL_GSCS),
                    UserWarning,
                    stacklevel=2,
                )
                self._is_custom = True
                logger.info(f"There is a custom GSC in {gsc_res!r}")
        self._gsc_res = gsc_res

    def load_genes(self, input_genes: List[str]):
        """Load gene list and convert to Entrez.

        Args:
            input_genes: Input gene list, can be mixed type.

        The following clsss attributes are set when ``load_genes`` is run

        :attr:`GenePlexus.input_genes` (List[str])
            Input genes converted to uppercase
        :attr:`GenePlexus.df_convert_out` (DataFrame)
            A table where the following 6 columns:

            .. list-table::

               * - Original ID
                 - User supplied Gene ID
               * - Entrez ID
                 - Entrez Gene ID
               * - Gene Name
                 - Name Gene ID
               * - In BioGRID?
                 - Y or N if the gene was found in the BioGRID network or not
               * - In IMP?
                 - Y or N if the gene was found in the IMP network or not
               * - In STRING?
                 - Y or N if the gene was found in the STRING network or not

        :attr:`GenePlexus.table_summary` (List[Dict[str, int]])
            List of netowrk stats summary dictionaries. Each dictionary has
            the following stucture:

            ::

               {
                 "Network" : # returns name of the network
                 "NetworkGenes"  : # returns number of genes in the network
                 "PositiveGenes" : # returns number of input genes found in the network
               }

        :attr:`GenePlexus.convert_ids` (List[str])
            Converted gene list.
        :attr:`GenePlexus.input_count` (int)
            Number of input genes that were able to be converted.

        See also:
            Use :meth:`geneplexus.util.read_gene_list` to load a gene list
            from a file.

        """
        self.input_genes = self._load_genes(input_genes)
        load_genes_outputs = self._convert_to_entrez(self.input_genes)
        self.df_convert_out = load_genes_outputs[0]
        self.table_summary = load_genes_outputs[1]
        self.input_count = load_genes_outputs[2]
        self.convert_ids = load_genes_outputs[3]
        self.model_info["All-Genes"].model_genes = self.convert_ids

    def load_negatives(self, input_negatives: List[str]):
        """Load gene list and convert to Entrez that will used as negatives.

        Args:
            input_negatives: Input negative gene list, can be mixed type.

        The following clsss attributes are set when ``load_negatives`` is run

        :attr:`GenePlexus.input_negatives` (List[str])
            Input negatives converted to uppercase
        :attr:`GenePlexus.df_convert_out_negatives` (DataFrame)
            A table with the following 6 columns:

            .. list-table::

               * - Original ID
                 - User supplied Gene ID
               * - Entrez ID
                 - Entrez Gene ID
               * - Gene Name
                 - Name Gene ID
               * - In BioGRID?
                 - Y or N if the gene was found in the BioGRID network or not
               * - In IMP?
                 - Y or N if the gene was found in the IMP network or not
               * - In STRING?
                 - Y or N if the gene was found in the STRING network or not

        :attr:`GenePlexus.table_summary_negatives` (List[Dict[str, int]])
            List of netowrk stats summary dictionaries. Each dictionary has
            the following stucture:

            ::

               {
                 "Network" : # returns name of the network
                 "NetworkGenes"  : # returns number of genes in the network
                 "PositiveGenes" : # returns number of input genes found in the network
               }

        :attr:`GenePlexus.convert_ids_negatives` (List[str])
            Converted negative gene list.
        :attr:`GenePlexus.input_count_negatives` (int)
            Number of negative genes that were able to be converted.

        See also:
            Use :meth:`geneplexus.util.read_gene_list` to load a gene list
            from a file.

        """
        self.input_negatives = self._load_genes(input_negatives)
        load_negatives_outputs = self._convert_to_entrez(self.input_negatives)
        self.df_convert_out_negatives = load_negatives_outputs[0]
        self.table_summary_negatives = load_negatives_outputs[1]
        self.input_count_negatives = load_negatives_outputs[2]
        self.convert_ids_negatives = load_negatives_outputs[3]

    def _load_genes(self, genes_to_load: List[str]):
        """Load gene list into the GenePlexus object.

        Note:
            Implicitely converts genes to upper case.

        """
        upper_genes = [item.upper() for item in genes_to_load]
        return upper_genes

    def _convert_to_entrez(self, genes_to_load: List[str]):
        """Convert the loaded genes to Entrez and make objects
        showing exactly what was converted

        """
        convert_ids, df_convert_out = _geneplexus._initial_id_convert(
            genes_to_load,
            self.file_loc,
            self.sp_trn,
        )
        df_convert_out, table_summary, input_count = _geneplexus._make_validation_df(
            df_convert_out,
            self.file_loc,
            self.sp_trn,
        )
        load_outputs = [df_convert_out, table_summary, input_count, convert_ids]
        return load_outputs

    def cluster_input(
        self,
        clust_method: config.CLUSTERING_TYPE = config.DEFAULT_PARAMETERS["clust_method"],
        clust_min_size: int = config.DEFAULT_PARAMETERS["clust_min_size"],
        clust_weighted: bool = config.DEFAULT_PARAMETERS["clust_weighted"],
        clust_kwargs: Optional[Dict[str, Any]] = config.DEFAULT_LOUVAIN_KWARGS | config.DEFAULT_DOMINO_KWARGS,
    ):
        """Cluster input gene list.

        Args:
            clust_method: Clustering method to use (either louvain or domino).
            clust_min_size: Ignore clusters if smaller than this value.
            clust_weighted: Whether or not to use weighted edges when building the clusters
            clust_kwargs: keywords args specfic to each clustering method
            louvain_max_size: (clust_kwarg, int) Try to recluster if a cluster is bigger than this value.
            louvain_max_tries: (clust_kwarg, int) The number of times to recluster any clusters that are
                bigger the `clust_max_size`. If cannot accomplished this by `clust_max_tries`
                the larger clusters are still retained.
            louvain_res: (clust_kwarg, float) Resolution parameter in clustering algorithm.
            louvain_seed: (clust_kwarg, int) Set seed used in clustering. Chose None to have this randomally set.
            domino_res: (clust_kwarg, float) resolution used to make initial slices.
            domino_slice_thresh: (clust_kwarg, float) threshold used for calling slice significant
            domino_n_steps: (clust_kwarg, int) number of steps used in pcst
            domino_module_threshold: (clust_kwarg, float) threshold used to consider module signifianct
            domino_seed: (clust_kwarg, int) random seed to be used in clustering algorithm

        The following clsss attributes are set when ``cluster_input`` is run

        :attr:`GenePlexus.clust_method` (str)
            Clustering method used
        :attr:`GenePlexus.clust_min_size` (int)
            Minimum size of clusters allowed
        :attr:`GenePlexus.clust_weighted` (bool)
            Whether or not to use edge weights when generating clusters
        :attr:`GenePlexus.clust_kwags` (dict)
            Keyword arguments used for each clustering method

        :attr:`GenePlexus.num_genes_lost` (int)
            Number of input_genes not in any cluster
        :attr:`GenePlexus.per_genes_lost` (float)
            Percentage of input_genes not in any cluster
        :attr:`GenePlexus.num_genes_gained` (int)
            Number of genes in clusters not in input_genes
        :attr:`GenePlexus.per_genes_gained` (float)
            Percentage of genes in clusters not in input_genes
        :attr:`GenePlexus.genes_lost_clustered` (List[str])
            List of input_genes not in any cluster
        :attr:`GenePlexus.genes_gained_clustered` (List[str])
            List of cluster genes not in input_genes
        :attr:`GenePlexus.model_info[ModelName].model_genes` (List[str])
            List of genes used as positives for each clusters model
        :attr:`GenePlexus.model_info[ModelName].results[ResultName]` (Class)
            For each clusters model, set up a key in results dicts for ModelResults class

        """

        if list(self.model_info) != ["All-Genes"]:
            warnings.warn(
                "\nOeleting previously generated clusters.",
                UserWarning,
                stacklevel=2,
            )
            for item in list(self.model_info):
                if "Cluster-" in item:
                    del self.model_info[item]

        if clust_method == "louvain":
            preset_kwargs = config.DEFAULT_LOUVAIN_KWARGS
        elif clust_method == "domino":
            preset_kwargs = config.DEFAULT_DOMINO_KWARGS
        preset_kwargs_keys = list(preset_kwargs.keys())
        if isinstance(clust_kwargs, dict):
            clust_kwargs = {key: value for key, value in clust_kwargs.items() if key in preset_kwargs_keys}
            preset_kwargs.update(clust_kwargs)
        else:
            logger.warning(f"clust_kwargs not a dictionary or None, using defaults")

        clust_genes = _geneplexus._generate_clusters(
            self.file_loc,
            self.sp_trn,
            self.net_type,
            self.model_info["All-Genes"].model_genes,
            clust_method,
            clust_min_size,
            clust_weighted,
            **preset_kwargs,
        )

        # set params to self for saving later
        self.clust_method = clust_method
        self.clust_min_size = clust_min_size
        self.clust_weighted = clust_weighted
        self.clus_kwargs = preset_kwargs
        # add keys to model_info
        if len(clust_genes) == 0:
            logger.info(f"No clusters were added")
        else:
            # add keys to model info
            clust_genes.sort(key=len, reverse=True)  # make biggest clusters first
            for i in range(len(clust_genes)):
                clus_id = i + 1
                self.model_info[f"Cluster-{clus_id:02d}"] = ModelInfo()
                self.model_info[f"Cluster-{clus_id:02d}"].model_genes = clust_genes[i]
                self.model_info[f"Cluster-{clus_id:02d}"].results = {}
                for apair in self.sp_gsc_pairs:
                    self.model_info[f"Cluster-{clus_id:02d}"].results[apair] = ModelResults()
            # generate info about the clusters
            unique_clus_genes = list({item for sublist in clust_genes for item in sublist})
            num_genes_lost = len(np.setdiff1d(self.model_info["All-Genes"].model_genes, unique_clus_genes))
            per_genes_lost = (num_genes_lost / len(self.model_info["All-Genes"].model_genes)) * 100
            num_genes_gained = len(np.setdiff1d(unique_clus_genes, self.model_info["All-Genes"].model_genes))
            per_genes_gained = (num_genes_gained / len(self.model_info["All-Genes"].model_genes)) * 100
            # set values for saving later
            self.num_genes_lost = num_genes_lost
            self.per_genes_lost = per_genes_lost
            self.num_genes_gained = num_genes_gained
            self.per_genes_gained = per_genes_gained
            self.genes_lost_clustered = np.setdiff1d(
                self.model_info["All-Genes"].model_genes,
                unique_clus_genes,
            ).tolist()
            self.genes_gained_clustered = np.setdiff1d(
                unique_clus_genes,
                self.model_info["All-Genes"].model_genes,
            ).tolist()
            logger.info(
                f"The number of clusters added is {len(clust_genes)}.\n"
                f"The number(%) of input genes removed by clustering is {num_genes_lost} ({per_genes_lost:.2f}%)\n"
                f"The number(%) of non-input genes added as positives by clustering is {num_genes_gained} ({per_genes_gained:.2f}%)",
            )

    def fit(
        self,
        logreg_kwargs: Optional[Dict[str, Any]] = config.DEFAULT_LOGREG_KWARGS,
        scale: bool = config.DEFAULT_PARAMETERS["scale"],
        min_num_pos: int = config.DEFAULT_PARAMETERS["min_num_pos"],
        min_num_pos_cv: int = config.DEFAULT_PARAMETERS["min_num_pos_cv"],
        num_folds: int = config.DEFAULT_PARAMETERS["num_folds"],
        null_val: float = config.DEFAULT_PARAMETERS["null_val"],
        random_state: Optional[int] = config.DEFAULT_PARAMETERS["random_state"],
        cross_validate: bool = config.DEFAULT_PARAMETERS["cross_validate"],
    ):
        """Fit the model.

        Args:
            logreg_kwargs: Scikit-learn logistic regression settings (see
                :class:`~sklearn.linear_model.LogisticRegression`).
            scale: Whether to scale the data when doing model training and prediction. It is
                not recommended to set to ``True`` unless using custom data.
            min_num_pos: Minimum number of positives required for the model
                to be trained.
            min_num_pos_cv: Minimum number of positives required for performing
                cross validation evaluation.
            num_folds: Number of cross validation folds.
            null_val: Null values to fill if cross validation was not able to
                be performed.
            random_state: Random state for reproducible shuffling stratified
                cross validation. Set to None for random.
            cross_validate: Whether or not to perform cross validation to
                evaluate the prediction performance on the gene set. If set to
                ``False``, then skip cross validation and return null_val as cv
                scores.

        The following clsss attributes are set when ``fit`` is run

        :attr:`GenePlexus.min_num_pos` (int)
            Minumum number of postivies needed to train a model.
        :attr:`GenePlexus.logreg_kwargs` (dict)
            Keyword arguments for LogisitcRegression function.
        :attr:`GenePlexus.scale` (bool)
            Whether or not scaling of the data was done in LogisticRegression.
        :attr:`GenePlexus.min_num_pos_cv` (int)
            The minumum number of positive genes needed for doing cross validation.
        :attr:`GenePlexus.num_folds` (int)
            Number of cross validation folds to do
        :attr:`GenePlexus.null_vall` (None, str, int, float)
            Value to fill in for avgps if cross validation couldn't be performed
        :attr:`GenePlexus.random_state` (None, int)
            Seed set for doing cross validation
        :attr:`GenePlexus.cross_validate` (bool)
            Whether or not to perform cross validation
        :attr:`GenePlexus.model_info[ModelName].pos_genes_in_net` (1D array of str)
            Input gene Entrez IDs that are present in the network.
        :attr:`GenePlexus.model_info[ModelName].genes_not_in_net` (1D array of str)
            Input gene Entrez IDs that are absent in the network.
        :attr:`GenePlexus.model_info[ModelName].net_genes` (1D array of str)
            All genes in the network.
        :attr:`GenePlexus.model_info[ModelName].negative_genes` (1D array of str)
            Negative gene Entrez IDs derived using the input genes and
            the background gene set collection (gp_trn).
        :attr:`GenePlexus.model_info[ModelName].neutral_gene_info` (Dict of Dicts)
            Dictionary saying which genes were set to neutrals because the
            term annotation matched closely enough to the positive training genes.

            ::

               {
                 "{Term ID}" # ID of the matched term : {
                    "Name"  : # returns string of term name
                    "Genes" : # returns list of genes annotated to term
                    "Task"  : # returns type of GSC the term is from
                    }
                 "All Neutrals" : # returns list of all genes considered neutral
               }

        :attr:`GenePlexus.model_info[ModelName].mdl_weights` (1D array of floats)
            Trained model parameters.
        :attr:`GenePlexus.model_info[ModelName].clf` (LogisticRegression)
            The fit classifer from sci-kit learn LogisticRegression class.
        :attr:`GenePlexusmodel_info[ModelName]..avgps` (1D array of floats)
            Cross validation results. Performance is measured using
            log2(auprc/prior).
        :attr:`GenePlexus.model_info[ModelName].std_scale` (StandardScale)
            If scaling was performed the object returned from StandardScaler.
        :attr:`GenePlexus.model_info[ModelName].df_convert_out_for_model` (DataFrame)
            A table specifc to input_genes for each model with the following 4 columns:

            .. list-table::

               * - Original ID
                 - User supplied Gene ID used to train the model
               * - Entrez ID
                 - Entrez Gene ID
               * - Gene Name
                 - Name Gene ID
               * - In {Network}?
                 - Y or N if the gene was found in the {Network} used to train the model

        Note:
            If setting scale to ``True`` then comparison of user trained model
            to the models pre-trained on known gene sets become less straightforward
            as those models are trained without any scaling.

        """
        if self.input_genes == None:
            raise NoPositivesError(
                f"No positives genes were added, use function load_genes()",
            )
        self.min_num_pos = min_num_pos
        if (self.clust_min_size != None) and (self.clust_min_size > self.min_num_pos):
            self.min_num_pos = self.clust_min_size
            logger.warning(
                "Setting the minimum number of genes to train a model to match the minumum allowable cluster size.",
            )

        for model_name in list(self.model_info):
            logger.info(f"Starting model training for {model_name}")
            self._get_pos_and_neg_genes(model_name)
            (
                self.model_info[model_name].mdl_weights,
                self.model_info[model_name].avgps,
                self.model_info[model_name].clf,
                self.model_info[model_name].std_scale,
            ) = _geneplexus._run_sl(
                self.file_loc,
                self.sp_trn,
                self.net_type,
                self.features,
                self.model_info[model_name].pos_genes_in_net,
                self.model_info[model_name].negative_genes,
                self.model_info[model_name].net_genes,
                logreg_kwargs,
                min_num_pos_cv,
                num_folds,
                null_val,
                random_state,
                cross_validate,
                scale,
            )

            # make df for genes used in training
            self.model_info[model_name].df_convert_out_for_model = _geneplexus._alter_validation_df(
                self.df_convert_out,
                self.model_info[model_name].pos_genes_in_net,
                self.net_type,
            )

        # set function arguments for saving later
        self.logreg_kwargs = logreg_kwargs
        self.scale = scale
        self.min_num_pos_cv = min_num_pos_cv
        self.num_folds = num_folds
        self.null_val = null_val
        self.random_state = random_state
        self.cross_validate = cross_validate
        return self.model_info

    def _get_pos_and_neg_genes(self, model_name):
        """Set up positive and negative splits."""
        (
            self.model_info[model_name].pos_genes_in_net,
            self.model_info[model_name].genes_not_in_net,
            self.model_info[model_name].net_genes,
        ) = _geneplexus._get_genes_in_network(
            self.file_loc,
            self.sp_trn,
            self.net_type,
            self.model_info[model_name].model_genes,
        )
        if len(self.model_info[model_name].pos_genes_in_net) < self.min_num_pos:
            raise NoPositivesError(
                f"There were not enough positive genes to train the model {model_name} with. "
                f"This limit is set to {self.min_num_pos} and can be changed in fit().",
            )

        if (self.input_negatives == None) or (len(self.input_negatives) == 0):
            user_negatives = None
        else:
            # remove genes from negatives if they are also positives
            user_negatives = np.setdiff1d(self.convert_ids_negatives, self.model_info[model_name].model_genes).tolist()
        (
            self.model_info[model_name].negative_genes,
            self.model_info[model_name].neutral_gene_info,
        ) = _geneplexus._get_negatives(
            self.file_loc,
            self.sp_trn,
            self.net_type,
            self.gsc_trn,
            self.model_info[model_name].pos_genes_in_net,
            user_negatives,
        )

        return (
            self.model_info[model_name].pos_genes_in_net,
            self.model_info[model_name].negative_genes,
            self.model_info[model_name].net_genes,
            self.model_info[model_name].neutral_gene_info,
        )

    def predict(self):
        """Predict gene scores from fit model.

        The following clsss attributes are set when ``predict`` is run

        :attr:`GenePlexus.model_info[ModelName].results[ResultName].df_probs` (DataFrame)
            A table with the following 9 columns:

            .. list-table::

               * - Entrez
                 - Entrez Gene ID
               * - Symbol
                 - Symbol Gene ID
               * - Name
                 - Name Gene ID
               * - Known/Novel
                 - Known is gene was in the positive set, otherwise Novel
               * - Class-Label
                 - P (positive in training), N (negative durinig training), U (unused during trianing)
               * - Probability
                 - The probabilties returned from the logisitc regression model
               * - Z-score
                 - The z-score of the model probabilties for all predcited genes
               * - P-adjusted
                 - The Bonferroni adjusted p-values from the z-scores
               * - Rank
                 - The rank of the gene with one being the gene with the highest predcited value


        Note:
            For the Known/Novel and Class-Label columns, if the training species is
            different than the results species, this information is obtained by looking
            at the one-to-one orthologs between the species.

        Note:
            Due to the high complexity of the embedding space, and wide variety of
            postive and negative genes determined for each model, the resulting
            probabilities may not be well calibrated, however the resulting rankings
            are very meaningful as evaluated with log2(auPRC/prior).

        """
        for model_name in list(self.model_info):
            for res_combo in list(self.model_info[model_name].results):
                logger.info(f"Generating predictions for {model_name} and {res_combo}")
                probs = _geneplexus._get_predictions(
                    self.file_loc,
                    res_combo.split("-")[0],
                    self.features,
                    self.net_type,
                    self.scale,
                    self.model_info[model_name].std_scale,
                    self.model_info[model_name].clf,
                )
                df_probs = _geneplexus._make_prob_df(
                    self.file_loc,
                    self.sp_trn,
                    res_combo.split("-")[0],
                    self.net_type,
                    probs,
                    self.model_info[model_name].pos_genes_in_net,
                    self.model_info[model_name].negative_genes,
                )
                self.model_info[model_name].results[res_combo].df_probs = df_probs
        return self.model_info

    def make_sim_dfs(self):
        """Compute similarities bewteen the input genes and GO, Monarch and/or Mondo.

        The following clsss attributes are set when ``make_sim_df`` is run

        :attr:`GenePlexus.model_info[ModelName].results[ResultName].df_sim` (DataFrame)
            A table showing how similar the coefficients of the user trained models
            are to the coefficients of models trained using genes annotated to gsc_res.
            The table has the following 7 columns:

            .. list-table::

               * - Task
                 - Which type of GSC the term is from
               * - ID
                 - Term ID
               * - Name
                 - Term Name
               * - Similarity
                 - Cosine similarity between model coefficients between the two models
               * - Z-score
                 - The z-score of the similarities
               * - P-adjusted
                 - The Bonferroni adjusted p-values from the z-scores
               * - Rank
                 - The rank of the term with one being the term with the highest similarity to the user model

        """
        for model_name in list(self.model_info):
            for idx, res_combo in enumerate(list(self.model_info[model_name].results)):
                logger.info(f"Generating model similarities for {model_name} and {res_combo}")
                df_sim, weights_dict = _geneplexus._make_sim_dfs(
                    self.file_loc,
                    self.model_info[model_name].mdl_weights,
                    res_combo.split("-")[0],
                    self.gsc_res[idx],  # needs to be different for Combines becoming GOs
                    self.net_type,
                    self.features,
                )
                self.model_info[model_name].results[res_combo].df_sim = df_sim
        return self.model_info

    def make_small_edgelist(
        self,
        num_nodes: int = config.DEFAULT_PARAMETERS["num_nodes"],
    ):
        """Make a subgraph induced by the top predicted genes.

        Args:
            num_nodes: Number of top genes to include.

        The following clsss attributes are set when ``make_small_edgelist`` is run

        :attr:`GenePlexus.num_nodes` (int)
            The number of nodes to include in the edgelist.
        :attr:`GenePlexus.model_info[ModelName].results[ResultName].df_edge` (DataFrame)
            Table of edge list corresponding to the subgraph induced by the top
            predicted genes (in Entrez gene ID).
        :attr:`GenePlexus.model_info[ModelName].results[ResultName].isolated_genes` (List[str])
            List of top predicted genes (in Entrez gene ID) that are isolated
            from other top predicted genes in the network.
        :attr:`GenePlexus.model_info[ModelName].results[ResultName].df_edge_sym` (DataFrame)
            Table of edge list corresponding to the subgraph induced by the top
            predicted genes (in gene symbol).
        :attr:`GenePlexus.model_info[ModelName].results[ResultName].isolated_genes_sym` (List[str])
            List of top predicted genes (in gene symbol) that are isolated from
            other top predicted genes in the network.

        """
        for model_name in list(self.model_info):
            for res_combo in list(self.model_info[model_name].results):
                logger.info(f"Generating small edgelists for {model_name} and {res_combo}")
                df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = _geneplexus._make_small_edgelist(
                    self.file_loc,
                    self.model_info[model_name].results[res_combo].df_probs,
                    res_combo.split("-")[0],
                    self.net_type,
                    num_nodes=num_nodes,
                )
                self.model_info[model_name].results[res_combo].df_edge = df_edge
                self.model_info[model_name].results[res_combo].isolated_genes = isolated_genes
                self.model_info[model_name].results[res_combo].df_edge_sym = df_edge_sym
                self.model_info[model_name].results[res_combo].isolated_genes_sym = isolated_genes_sym
        # set value for saving later
        self.num_nodes = num_nodes
        return self.model_info

    def save_class(
        self,
        output_dir: str = config.DEFAULT_PARAMETERS["output_dir"],
        save_type: config.SAVE_TYPE = config.DEFAULT_PARAMETERS["save_type"],
        zip_output: bool = config.DEFAULT_PARAMETERS["zip_output"],
        overwrite: bool = config.DEFAULT_PARAMETERS["overwrite"],
    ):
        """Save all or parts of the GenePlexus class and results.

        Args:
            output_dir: Path to save the files to If None will try ~/.data/geneplexus_outputs/results.
            save_type: which file saving method to use
            zip_output: wehter or not to compress all the results into one zip file
            overwrite: wether to overwrite data or make new directory with incremented index

        """

        _geneplexus._save_class(self, output_dir, save_type, zip_output, overwrite)

    def remove_log_file(self):
        """Remove the tmp log file. Only do when at the end of the script)"""

        if self.log_to_file:
            if os.path.exists(self.log_tmp_path):
                logger.removeHandler(self.file_handler)
                os.remove(self.log_tmp_path)
