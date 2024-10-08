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
from .exception import CustomDataError
from .exception import FlyMonarchError
from .exception import ZebrafishBioGRIDError


class GenePlexus:
    """The GenePlexus API class."""

    def __init__(
        self,
        file_loc: Optional[str] = None,
        net_type: config.NET_TYPE = "STRING",
        features: config.FEATURE_TYPE = "SixSpeciesN2V",
        sp_trn: config.SPECIES_TYPE = "Human",
        sp_tst: config.SPECIES_TYPE = "Human",
        gsc_trn: config.GSC_TYPE = "Combined",
        gsc_tst: config.GSC_TYPE = "Combined",
        input_genes: Optional[List[str]] = None,
        input_negatives: Optional[List[str]] = None,
        auto_download: bool = False,
        log_level: config.LOG_LEVEL_TYPE = "WARNING",
    ):
        """Initialize the GenePlexus object.

        Args:
            file_loc: Location of data files, if not specified, set to default
                data path ``~/.data/geneplexus``
            net_type: Type of network to use.
            features: Type of features of the network to use.
            sp_trn: The species of the training data
            sp_tst: The species of the testing data
            gsc_trn: Gene set collection used during training
            gsc_tst: Gene set collection used when generating
                similarity dataframe
            input_genes: Input gene list, can be mixed type. Can also be set
                later if not specified at init time by simply calling
                :meth:`load_genes` (default: :obj:`None`).
            input_negatives: Input list of negative genes, can be mixed type.
                Can also be set later if not specified at init time by simply calling
                :meth:`load_negatives` (default: :obj:`None`).
            auto_download: Automatically download necessary files if set.
            log_level: Logging level.

        """
        set_stream_level(logger, log_level)
        self._is_custom: bool = False
        self.file_loc = file_loc  # type: ignore
        self.features = features
        self.sp_trn = sp_trn
        self.sp_tst = sp_tst
        self.gsc_trn = gsc_trn
        self.gsc_tst = gsc_tst
        self.net_type = net_type
        self.log_level = log_level
        self.auto_download = auto_download
        self.input_genes: List[str] = []
        self.input_negatives: List[str] = []

        if self.auto_download and self._is_custom:
            warnings.warn(
                "\nSkipping auto download for custom files. Unset auto_download option to suppress this message.",
                UserWarning,
                stacklevel=2,
            )
        elif self.auto_download:
            download_select_data(
                self.file_loc,
                list({self.sp_trn, self.sp_tst}),
                log_level=log_level,
            )

        if input_genes is not None:
            self.load_genes(input_genes)

        if input_negatives is not None:
            self.load_negatives(input_negatives)

        if ("Zebrafish" == (self.sp_trn or self.sp_tst)) and (self.net_type == "BioGRID"):
            raise ZebrafishBioGRIDError(
                f"The BioGRID network for Zebrafish is not "
                "included due to it not having enough nodes "
                "so this combination is not allowed.",
            )

        if (
            (self.sp_trn == "Fly" and self.gsc_trn == "Monarch")
            or (self.sp_trn == "Fly" and self.gsc_trn == "Combined")
            or (self.sp_tst == "Fly" and self.gsc_tst == "Monarch")
            or (self.sp_tst == "Fly" and self.gsc_tst == "Combined")
        ):
            raise FlyMonarchError(
                f"Fly has no annotations for Monarch. Use GO for GSC",
            )

    @property
    def _params(self) -> List[str]:
        return [
            "file_loc",
            "net_type",
            "features",
            "sp_trn",
            "sp_tst",
            "gsc_trn",
            "gsc_tst",
            "auto_download",
            "log_level",
            "input_genes",
            "input_negatives",
        ]

    def dump_config(self, outdir: str):
        """Save parameters configuration to a config file
        when running with CLI.
        """
        params_dict = {i: getattr(self, i) for i in self._params}
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
        """Geneset collection."""
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
    def sp_tst(self) -> config.SPECIES_TYPE:
        """Geneset collection."""
        return self._sp_tst

    @sp_tst.setter
    def sp_tst(self, sp_tst: config.SPECIES_TYPE):
        if sp_tst not in config.ALL_SPECIES:
            warnings.warn(
                util.param_warning("species", sp_tst, config.ALL_SPECIES),
                UserWarning,
                stacklevel=2,
            )
            self._is_custom = True
            logger.info(f"Using custom species {sp_tst!r}")
        self._sp_tst = sp_tst

    @property
    def gsc_trn(self) -> config.GSC_TYPE:
        """Geneset collection."""
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
    def gsc_tst(self) -> config.GSC_TYPE:
        """Geneset collection."""
        return self._gsc_tst

    @gsc_tst.setter
    def gsc_tst(self, gsc_tst: config.GSC_TYPE):
        if gsc_tst not in config.ALL_GSCS:
            warnings.warn(
                util.param_warning("GSC", gsc_tst, config.ALL_GSCS),
                UserWarning,
                stacklevel=2,
            )
            self._is_custom = True
            logger.info(f"Using custom GSC {gsc_tst!r}")
        self._gsc_tst = gsc_tst

    def load_genes(self, input_genes: List[str]):
        """Load gene list and convert to Entrez.

        :attr:`GenePlexus.input_genes` (List[str]): Input gene list.

        Args:
            input_genes: Input gene list, can be mixed type.

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

    def load_negatives(self, input_negatives: List[str]):
        """Load gene list and convert to Entrez that will used as negatives.

        :attr:`GenePlexus.input_negatives` (List[str]): Input gene list.

        Args:
            input_negstives: Negative gene list, can be mixed type.

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
        """Convert the loaded genes to Entrez.

        :attr:`GenePlexus.df_convert_out` (DataFrame)
            A table where the first column contains the original gene IDs, the
            second column contains the corresponding converted Entrez gene IDs.
            The rest of the columns are indicators of whether a given gene is
            present in any one of the networks.
        :attr:`GenePlexus.table_summary` (List[Dict[str, int]])
            List of netowrk stats summary dictionaries. Each dictionary has
            three keys: **Network**, **NetworkGenes**, and **PositiveGenes**
            (the number intersection between the input genes and the network
            genes).
        :attr:`GenePlexus.input_count` (int)
            Number of input genes.

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

    def fit_and_predict(
        self,
        logreg_kwargs: Optional[Dict[str, Any]] = None,
        min_num_pos: int = 15,
        num_folds: int = 3,
        null_val: float = None,
        random_state: Optional[int] = 0,
        cross_validate: bool = True,
    ):
        """Fit a model and predict gene scores.

        Args:
            logreg_kwargs: Scikit-learn logistic regression settings (see
                :class:`~sklearn.linear_model.LogisticRegression`). If not set,
                then use the default logistic regression settings (l2 penalty,
                10,000 max iterations, lbfgs solver).
            min_num_pos: Minimum number of positives required for performing
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

        :attr:`GenePlexus.mdl_weights` (array of float)
            Trained model parameters.
        :attr:`GenePlexus.probs` (array of float)
            Genome-wide gene prediction scores. A high value indicates the
            relevance of the gene to the input gene list.
        :attr:`GenePlexus.avgps` (array of float)
            Cross validation results. Performance is measured using
            log2(auprc/prior).
        :attr:`GenePlexus.df_probs` (DataFrame)
            A table with 7 columns: **Entrez** (the gene Entrez ID), **Symbol**
            (the gene Symbol), **Name** (the gene Name), **Probability** (the
            probability of a gene being part of the input gene list),
            **Known/Novel** (whether the gene is in the input gene list),
            **Class-Label** (positive, negative, or neutral), **Rank** (rank of
            relevance of the gene to the input gene list).

        """
        self._get_pos_and_neg_genes()
        self.mdl_weights, self.probs, self.avgps = _geneplexus._run_sl(
            self.file_loc,
            self.sp_trn,
            self.sp_tst,
            self.net_type,
            self.features,
            self.pos_genes_in_net,
            self.negative_genes,
            self.net_genes,
            logreg_kwargs=logreg_kwargs,
            min_num_pos=min_num_pos,
            num_folds=num_folds,
            null_val=null_val,
            random_state=random_state,
            cross_validate=cross_validate,
        )
        self.df_probs = _geneplexus._make_prob_df(
            self.file_loc,
            self.sp_trn,
            self.sp_tst,
            self.net_type,
            self.probs,
            self.pos_genes_in_net,
            self.negative_genes,
        )
        return self.mdl_weights, self.df_probs, self.avgps

    def _get_pos_and_neg_genes(self):
        """Set up positive and negative genes given the network.

        :attr:`GenePlexus.pos_genes_in_net` (array of str)
            Array of input gene Entrez IDs that are present in the network.
        :attr:`GenePlexus.genes_not_in_net` (array of str)
            Array of input gene Entrez IDs that are absent in the network.
        :attr:`GenePlexus.net_genes` (array of str)
            Array of network gene Entrez IDs.
        :attr:`GenePlexus.negative_genes` (array of str)
            Array of negative gene Entrez IDs derived using the input genes and
            the background gene set collection (GSC).

        """
        self.pos_genes_in_net, self.genes_not_in_net, self.net_genes = _geneplexus._get_genes_in_network(
            self.file_loc,
            self.sp_trn,
            self.net_type,
            self.convert_ids,
        )
        if len(self.input_negatives) > 0:
            # remove genes from negatives if they are also positives
            user_negatives = np.setdiff1d(self.convert_ids_negatives, self.convert_ids).tolist()
            print(user_negatives)
        else:
            user_negatives = None
        self.negative_genes, self.neutral_gene_info = _geneplexus._get_negatives(
            self.file_loc,
            self.sp_trn,
            self.net_type,
            self.gsc_trn,
            self.pos_genes_in_net,
            user_negatives,
        )

        return self.pos_genes_in_net, self.negative_genes, self.net_genes, self.neutral_gene_info

    def make_sim_dfs(self):
        """Compute similarities bewteen the input genes and GO or Mondo.

        The similarities are compuared based on the model trained on the input
        gene set and models pre-trained on known GO and Mondo gene sets.

        :attr:`GenePlexus.df_sim_GO` (DataFrame)
            A table with 4 columns: **ID** (the GO term ID), **Name** (name of
            the GO term), **Similarity** (similarity between the input model
            and a model trained on the GO term gene set), **Rank** (rank of
            similarity between the input model and a model trained on the GO
            term gene set).
        :attr:`GenePlexus.df_sim_Dis` (DataFrame)
            A table with 4 columns: **ID** (the DO term ID), **Name** (name of
            the DO term), **Similarity** (similarity between the input model
            and a model trained on the DO term gene set), **Rank** (rank of
            similarity between the input model and a model trained on the DO
            term gene set).
        :attr:`GenePlexus.weights_GO`
            Dictionary of pretrained model weights for GO. A key is a GO term,
            and the value is a dictionary with three keys: **Name** (name of
            the GO term), **Weights** (pretrained model weights), **PosGenes**
            (positive genes for this GO term).
        :attr:`GenePlexus.weights_Dis`
            Dictionary of pretrained model weights for Mondo. A key is a DO
            term, and the value is a dictionary with three keys: **Name** (name
            of the DO term), **Weights** (pretrained model weights),
            **PosGenes** (positive genes for this DO term).

        """
        self.df_sim, self.weights_dict = _geneplexus._make_sim_dfs(
            self.file_loc,
            self.mdl_weights,
            self.sp_tst,
            self.gsc_tst,
            self.net_type,
            self.features,
        )
        return self.df_sim, self.weights_dict

    def make_small_edgelist(self, num_nodes: int = 50):
        """Make a subgraph induced by the top predicted genes.

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

        Args:
            num_nodes: Number of top genes to include.

        """
        self.df_edge, self.isolated_genes, self.df_edge_sym, self.isolated_genes_sym = _geneplexus._make_small_edgelist(
            self.file_loc,
            self.df_probs,
            self.sp_tst,
            self.net_type,
            num_nodes=num_nodes,
        )
        return self.df_edge, self.isolated_genes, self.df_edge_sym, self.isolated_genes_sym

    def alter_validation_df(self):
        """Make table about presence of input genes in the network.

        :attr:`df_convert_out_subset`
        :attr:`positive_genes`

        """
        self.df_convert_out_subset, self.positive_genes = _geneplexus._alter_validation_df(
            self.df_convert_out,
            self.table_summary,
            self.net_type,
        )
        return self.df_convert_out_subset, self.positive_genes
