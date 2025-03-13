"""GenePlexus API."""
import os
import os.path as osp
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
        file_loc: Optional[str] = None,
        net_type: config.NET_TYPE = "STRING",
        features: config.FEATURE_TYPE = "SixSpeciesN2V",
        sp_trn: config.SPECIES_TYPE = "Human",
        sp_res: config.SPECIES_SELECTION_TYPE = "Human",
        gsc_trn: config.GSC_TYPE = "Combined",
        gsc_res: config.GSC_SELECTION_TYPE = "Combined",
        input_genes: Optional[List[str]] = None,
        input_negatives: Optional[List[str]] = None,
        auto_download: bool = False,
        log_level: config.LOG_LEVEL_TYPE = "WARNING",
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

        :attr:`GenePlexus.gsc_trn_original` str
            If internal data checks are run, this can different that gsc_trn.
        :attr:`GenePlexus.gsc_res_original` (List[str])
            If internal data checks are run, this can different that res_trn.

        """
        set_stream_level(logger, log_level)
        self._is_custom: bool = False
        self.file_loc = file_loc  # type: ignore
        self.features = features
        self.sp_trn = sp_trn
        self.sp_res = sp_res
        self.gsc_trn = gsc_trn
        self.gsc_res = gsc_res
        self.net_type = net_type
        self.log_level = log_level
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

        if input_genes is not None:
            self.load_genes(input_genes)

        if input_negatives is not None:
            self.load_negatives(input_negatives)

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

        # remove duplicate sp-gsc comboms if any for results
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

    @property
    def _params(self) -> List[str]:
        return [
            "file_loc",
            "net_type",
            "features",
            "sp_trn",
            "sp_res",
            "gsc_trn",
            "gsc_res",
            "auto_download",
            "log_level",
            "input_genes",
            "input_negatives",
        ]

    def dump_config(self, outdir: str):
        """Save parameters configuration to a config file,
        used with CLI.
        """
        params_dict = {i: getattr(self, i) for i in self._params}
        if params_dict["gsc_trn"] == "Combined":
            params_dict["gsc_trn"] = config.COMBINED_CONTEXTS[params_dict["sp_trn"]]
        if params_dict["gsc_res"] == "Combined":
            params_dict["gsc_res"] = config.COMBINED_CONTEXTS[params_dict["sp_res"]]
        if hasattr(self, "gsc_trn_original") and (self.gsc_trn_original == "Combined"):
            params_dict["gsc_trn"] = config.COMBINED_CONTEXTS[params_dict["sp_trn"]]
        if hasattr(self, "gsc_res_original") and (self.gsc_res_original == "Combined"):
            params_dict["gsc_res"] = config.COMBINED_CONTEXTS[params_dict["sp_res"]]
        path = osp.join(outdir, "config.yaml")
        with open(path, "w") as f:
            yaml.dump(params_dict, f)
            logger.info(f"Config saved to {path}")

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

        **The following clsss attributes are set when this function is run**

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

        **The following clsss attributes are set when this function is run**

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
        clust_min_size: int = 5,
        clust_max_size: int = 70,
        clust_max_tries: int = 3,
        clust_res: int = 1,
        clust_weighted: bool = True,
    ):
        """Cluster input gene list.

        Args:
            clust_min_size: Ignore clusters if smaller than this value.
            clust_max_size: Try to recluster if a cluster is bigger than this value.
            clust_max_tries: The number of times to recluster any clusters that are
                bigger the `clust_max_size`. If cannot accomplished this by `clust_max_tries`
                the larger clusters are still retained.
            clust_res: Resolution parameter in clustering algorithm.
            clust_weighted: Wether or not to use weighted edges when building the clusters
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

        clust_genes = _geneplexus._generate_clusters(
            self.file_loc,
            self.sp_trn,
            self.net_type,
            self.model_info["All-Genes"].model_genes,
            clust_min_size,
            clust_max_size,
            clust_max_tries,
            clust_res,
            clust_weighted,
        )
        # add keys to model_info
        if len(clust_genes) == 0:
            logger.info(f"No clusters were added")
        else:
            clust_genes.sort(key=len, reverse=True)  # make biggest clusters first
            for i in range(len(clust_genes)):
                clus_id = i + 1
                self.model_info[f"Cluster-{clus_id:02d}"] = ModelInfo()
                self.model_info[f"Cluster-{clus_id:02d}"].model_genes = clust_genes[i]
                self.model_info[f"Cluster-{clus_id:02d}"].results = {}
                for apair in self.sp_gsc_pairs:
                    self.model_info[f"Cluster-{clus_id:02d}"].results[apair] = ModelResults()

    def fit(
        self,
        logreg_kwargs: Optional[Dict[str, Any]] = None,
        scale: bool = False,
        min_num_pos: int = 5,
        min_num_pos_cv: int = 15,
        num_folds: int = 3,
        null_val: float = None,
        random_state: Optional[int] = 0,
        cross_validate: bool = True,
    ):
        """Fit the model.

        Args:
            logreg_kwargs: Scikit-learn logistic regression settings (see
                :class:`~sklearn.linear_model.LogisticRegression`). If not set,
                then use the default logistic regression settings (l2 penalty,
                10,000 max iterations, lbfgs solver).
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

        **The following clsss attributes are set when this function is run**

        :attr:`GenePlexus.mdl_weights` (1D array of floats)
            Trained model parameters.
        :attr:`GenePlexus.avgps` (1D array of floats)
            Cross validation results. Performance is measured using
            log2(auprc/prior).

        Note:
            If setting scale to ``True`` then comparison of user trained model
            to the models pre-trained on known gene sets become less straightforward
            as those models are trained without any scaling.

        """
        if self.input_genes == None:
            raise NoPositivesError(
                f"No positives genes were added, use function load_genes()",
            )
        for model_name in list(self.model_info):
            self._get_pos_and_neg_genes(model_name, min_num_pos)
            (
                self.model_info[model_name].mdl_weights,
                self.model_info[model_name].avgps,
                self.model_info[model_name].scale,
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
                logreg_kwargs=logreg_kwargs,
                min_num_pos_cv=min_num_pos_cv,
                num_folds=num_folds,
                null_val=null_val,
                random_state=random_state,
                cross_validate=cross_validate,
                scale=scale,
            )
        return self.model_info

    def _get_pos_and_neg_genes(self, model_name, min_num_pos):
        """Set up positive and negative splits.

        **The following clsss attributes are set when this function is run**

        :attr:`GenePlexus.pos_genes_in_net` (1D array of str)
            Input gene Entrez IDs that are present in the network.
        :attr:`GenePlexus.genes_not_in_net` (1D array of str)
            Input gene Entrez IDs that are absent in the network.
        :attr:`GenePlexus.net_genes` (1D array of str)
            All genes in the network.
        :attr:`GenePlexus.negative_genes` (1D array of str)
            Negative gene Entrez IDs derived using the input genes and
            the background gene set collection (gp_trn).
        :attr:`GenePlexus.neutral_gene_info` (Dict of Dicts)
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

        """
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
        if len(self.model_info[model_name].pos_genes_in_net) < min_num_pos:
            raise NoPositivesError(
                f"There were not enough positive genes to train the model {model_name} with. "
                f"This limit is set to {min_num_pos} and can be changed in fit().",
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

        **The following clsss attributes are set when this function is run**

        :attr:`GenePlexus.df_probs` (DataFrame)
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

        :attr:`GenePlexus.probs` (1D array of floats)
            Genome-wide gene prediction scores. A high value indicates the
            relevance of the gene to the input gene list.

        """
        for model_name in list(self.model_info):
            for res_combo in list(self.model_info[model_name].results):
                probs = _geneplexus._get_predictions(
                    self.file_loc,
                    res_combo.split("-")[0],
                    self.features,
                    self.net_type,
                    self.model_info[model_name].scale,
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
                self.model_info[model_name].results[res_combo].probs = probs
                self.model_info[model_name].results[res_combo].df_probs = df_probs
        return self.model_info

    def make_sim_dfs(self):
        """Compute similarities bewteen the input genes and GO, Monarch and/or Mondo.

        **The following clsss attributes are set when this function is run**

        :attr:`GenePlexus.df_sim` (DataFrame)
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

        :attr:`GenePlexus.weights`
            Dictionary of pretrained model weights for gsc_res.

            ::

               {
                 "{Term ID}" # ID of the GSC term : {
                    "Name"  : # returns string of term name
                    "PosGenes" : # returns list of genes annotated to term
                    "Task"  : # returns type of GSC the term is from
                    "Weights" : # return list of coefficients from models trained using genes annotated to the term
                    }
               }

        """
        for model_name in list(self.model_info):
            for idx, res_combo in enumerate(list(self.model_info[model_name].results)):
                df_sim, weights_dict = _geneplexus._make_sim_dfs(
                    self.file_loc,
                    self.model_info[model_name].mdl_weights,
                    res_combo.split("-")[0],
                    self.gsc_res[idx],  # needs to be different for Combines becoming GOs
                    self.net_type,
                    self.features,
                )
                self.model_info[model_name].results[res_combo].df_sim = df_sim
                self.model_info[model_name].results[res_combo].weights_dict = weights_dict
        return self.model_info

    def make_small_edgelist(self, num_nodes: int = 50):
        """Make a subgraph induced by the top predicted genes.

        Args:
            num_nodes: Number of top genes to include.

        **The following clsss attributes are set when this function is run**

        :attr:`GenePlexus.df_edge` (DataFrame)
            Table of edge list corresponding to the subgraph induced by the top
            predicted genes (in Entrez gene ID).
        :attr:`GenePlexus.isolated_genes` (List[str])
            List of top predicted genes (in Entrez gene ID) that are isolated
            from other top predicted genes in the network.
        :attr:`GenePlexus.df_edge_sym` (DataFrame)
            Table of edge list corresponding to the subgraph induced by the top
            predicted genes (in gene symbol).
        :attr:`GenePlexus.isolated_genes_sym` (List[str])
            List of top predicted genes (in gene symbol) that are isolated from
            other top predicted genes in the network.

        """
        for model_name in list(self.model_info):
            for res_combo in list(self.model_info[model_name].results):
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
        return self.model_info

    def alter_validation_df(self):
        """Make table about presence of input genes in the network used durning training.

        **The following clsss attributes are set when this function is run**

        :attr:`df_convert_out_subset` (DataFrame)
            A table with the following 6 columns:

            .. list-table::

               * - Original ID
                 - User supplied Gene ID
               * - Entrez ID
                 - Entrez Gene ID
               * - Gene Name
                 - Name Gene ID
               * - In {Network}?
                 - Y or N if the gene was found in the {Network} used to train the model

        :attr:`positive_genes` (List[str])
            List of genes used as positives when training the model


        """
        self.df_convert_out_subset, self.positive_genes = _geneplexus._alter_validation_df(
            self.df_convert_out,
            self.table_summary,
            self.net_type,
        )
        return self.df_convert_out_subset, self.positive_genes
