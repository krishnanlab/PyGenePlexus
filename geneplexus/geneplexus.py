from typing import List

from . import _geneplexus
from ._config import config


class GenePlexus:
    def __init__(
        self,
        file_loc: str,
        net_type: config.NET_TYPE = "BioGRID",
        features: config.FEATURE_TYPE = "Embedding",
        GSC: config.GSC_TYPE = "GO",
    ):
        """Initialize the GenePlexus object.

        Args:
            file_loc (str): Location of data files.
            net_type (NET_TYPE): Type of network to use.
            features (FEATURE_TYPE): Type of features of the network to use.
            GSC (GSC_TYPE): Type of gene set collection to use for generating
                negatives.

        Attributes:
            input_genes (List[str]): Input gene list. Created by
                :meth:`GenePlexus.load_genes`
            df_convert_out (DataFrame): A table where the first column contains
                the original gene IDs, the second column contains the
                corresponding converted Entrez gene IDs. The rest of the
                columns are indicators of whether a given gene is present in
                any one of the networks. Created by
                :meth:`GenePlexus.convert_to_Entrez`
            table_summary (List[Dict[str, int]]): List of netowrk stats summary
                dictionaries. Each dictionary has three keys: **Network**,
                **NetworkGenes**, and **PositiveGenes** (the number
                intersection between the input genes and the network genes).
                Created by :meth:`GenePlexus.convert_to_Entrez`
            input_count (int): Number of input genes. Created by
                :meth:`GenePlexus.convert_to_Entrez`
            pos_genes_in_net (array of str): Array of input gene Entrez IDs
                that are present in the network. Created by
                :meth:`GenePlexus.get_pos_and_neg_genes`
            genes_not_in_net (array of str): Array of input gene Entrez IDs
                that are absent in the network. Created by
                :meth:`GenePlexus.get_pos_and_neg_genes`
            net_genes (array of str): Array of network gene Entrez IDs.
                Created by :meth:`GenePlexus.get_pos_and_neg_genes`
            negative_genes (array of str): Array of negative gene Entrez IDs
                derived using the input genes and the background gene set
                collection (GSC). Created by
                :meth:`GenePlexus.get_pos_and_neg_genes`
            mdl_weights (array of float): Trained model parameters. Created
                by :meth:`GenePlexus.fit_and_predict`
            probs (array of float): Genome wide gene prediction scores. A high
                value indicates the relevance of the gene to the input gene
                list. Created by :meth:`GenePlexus.fit_and_predict`
            avgps (array of float): Cross validation results. Performance
                is measured using log2(auprc/prior). Created by
                :meth:`GenePlexus.fit_and_predict`
            df_probs (DataFrame): A table with 7 columns: **Entrez** (the gene
                Entrez ID), **Symbol** (the gene Symbol), **Name** (the gene
                Name), **Probability** (the probability of a gene being part of
                the input gene list), **Known/Novel** (whether the gene is in
                the input gene list), **Class-Label** (positive, negative, or
                neutral), **Rank** (rank of relevance of the gene to the input
                gene list). Created by :meth:`GenePlexus.fit_and_predict`
            df_sim_GO (DataFrame): A table with 4 columns: **ID** (the GO term
                ID), **Name** (name of the GO temr), **Similarity** (similarity
                between the input gene list and a GO term), **Rank** (rank of
                GO term similarity with the input gene list). Created by
                :meth:`GenePlexus.make_sim_dfs`
            df_sim_Dis (DataFrame): A table with 4 columns: **ID** (the DO term
                ID), **Name** (name of the DO temr), **Similarity** (similarity
                between the input gene list and a DO term), **Rank** (rank of
                DO term similarity with the input gene list). Created by
                :meth:`GenePlexus.make_sim_dfs`
            weights_GO: Dictionary of pretrained model weights for GO. A key is
                a GO term, and the value is a dictionary with three keys:
                **Name** (name of the GO term), **Weights** (pretrained model
                weights), **PosGenes** (positive genes for this GO term).
                Created by :meth:`GenePlexus.make_sim_dfs`
            weights_Dis: Dictionary of pretrained model weights for DisGeNet. A
                key is a DO term, and the value is a dictionary with three
                keys: **Name** (name of the DO term), **Weights** (pretrained
                model weights), **PosGenes** (positive genes for this DO term).
                Created by :meth:`GenePlexus.make_sim_dfs`
            df_edge (DataFrame): Table of edge list corresponding to the
                subgraph induced by the top predicted genes (in Entrez gene
                ID). Created by :meth:`make_small_edgelist`
            isolated_genes (List[str]): List of top predicted genes (in Entrez
                gene ID) are are isolated from other top predicted genes in
                the network. Created by :meth:`make_small_edgelist`
            df_edge_sym (DataFrame): Table of edge list corresponding to the
                subgraph induced by the top predicted genes (in gene symbol).
                Created by :meth:`make_small_edgelist`
            isolated_genes_sym (List[str]): List of top predicted genes (in gene
                symbol) are are isolated from other top predicted genes in the
                network. Created by :meth:`make_small_edgelist`
            df_convert_out_subset: Created by :meth:`alter_validation_df`
            positive_genes: Created by :meth:`alter_validation_df`

        Todos:
            - :attr:`genes_not_in_net` has wrong type (array of int).

        """
        self.file_loc = file_loc
        self.net_type = net_type
        self.features = features
        self.GSC = GSC

    def load_genes(self, input_genes: List[str]):
        """Load list of genes into the GenePlexus object.

        Creates :attr:`input_genes`

        Args:
            input_genes (List[str]): Input genes, can be mixed type.

        See also:
            Use :meth:`geneplexus.util.read_gene_list` to load a gene list
            from a file.

        """
        input_genes = [item.upper() for item in input_genes]
        self.input_genes = input_genes

    def convert_to_Entrez(self):
        """Convert the loaded genes to Entrez.

        Creates :attr:`df_convert_out`, :attr:`table_summary`, and
        :attr:`input_count`

        """
        self.convert_IDs, df_convert_out = _geneplexus.initial_ID_convert(self.input_genes, self.file_loc)
        self.df_convert_out, self.table_summary, self.input_count = _geneplexus.make_validation_df(
            df_convert_out,
            self.file_loc,
        )
        return self.df_convert_out

    def set_params(
        self,
        net_type: config.NET_TYPE,
        features: config.FEATURE_TYPE,
        GSC: config.GSC_TYPE,
    ):
        """Set GenePlexus parameters."""
        self.net_type = net_type
        self.features = features
        self.GSC = GSC

    def get_pos_and_neg_genes(self):
        """Set up positive and negative genes given the network.

        Creates :attr:`pos_genes_in_net`, :attr:`genes_not_in_net`,
        :attr:`net_genes`, :attr:`negative_genes`

        """
        self.pos_genes_in_net, self.genes_not_in_net, self.net_genes = _geneplexus.get_genes_in_network(
            self.file_loc,
            self.net_type,
            self.convert_IDs,
        )
        self.negative_genes = _geneplexus.get_negatives(
            self.file_loc,
            self.net_type,
            self.GSC,
            self.pos_genes_in_net,
        )
        return self.pos_genes_in_net, self.negative_genes, self.net_genes

    def fit_and_predict(self):
        """Fit a model and predict gene scores.

        Creates :attr:`mdl_weights`, :attr:`probs`, :attr:`avgps`, and
        :attr:`df_probs`

        """
        self.mdl_weights, self.probs, self.avgps = _geneplexus.run_SL(
            self.file_loc,
            self.net_type,
            self.features,
            self.pos_genes_in_net,
            self.negative_genes,
            self.net_genes,
        )
        self.df_probs = _geneplexus.make_prob_df(
            self.file_loc,
            self.net_genes,
            self.probs,
            self.pos_genes_in_net,
            self.negative_genes,
        )
        return self.mdl_weights, self.df_probs, self.avgps

    def make_sim_dfs(self):
        """Compute similarities bewteen the input genes and GO or DisGeNet.

        Creates :attr:`df_sim_GO`, :attr:`df_sim_Dis`, :attr:`weights_GO`, and
        :attr:`weights_Dis`

        """
        self.df_sim_GO, self.df_sim_Dis, self.weights_GO, self.weights_Dis = _geneplexus.make_sim_dfs(
            self.file_loc,
            self.mdl_weights,
            self.GSC,
            self.net_type,
            self.features,
        )
        return self.df_sim_GO, self.df_sim_Dis, self.weights_GO, self.weights_Dis

    def make_small_edgelist(self, num_nodes: int = 50):
        """Make a subgraph induced by the top predicted genes.

        Creates :attr:`df_edge`, :attr:`isolated_genes`, :attr:`df_edge_sym`,
        and :attr:`isolated_genes_sym`

        Args:
            num_nodes (int): Number of top genes to include.

        """
        self.df_edge, self.isolated_genes, self.df_edge_sym, self.isolated_genes_sym = _geneplexus.make_small_edgelist(
            self.file_loc,
            self.df_probs,
            self.net_type,
            num_nodes=50,
        )
        return self.df_edge, self.isolated_genes, self.df_edge_sym, self.isolated_genes_sym

    def alter_validation_df(self):
        """Make table about presence of input genes in the network.

        Creates :attr:`df_convert_out_subset` and :attr:`positive_genes`

        """
        self.df_convert_out_subset, self.positive_genes = _geneplexus.alter_validation_df(
            self.df_convert_out,
            self.table_summary,
            self.net_type,
        )
        return self.df_convert_out_subset, self.positive_genes
