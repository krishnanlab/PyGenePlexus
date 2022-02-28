from . import _geneplexus


class GenePlexus:
    def __init__(
        self,
        file_loc,
        net_type="BioGRID",
        features="Embedding",
        GSC="GO",
    ):
        self.file_loc = file_loc
        self.net_type = net_type
        self.features = features
        self.GSC = GSC

    def load_genes(self, input_genes):
        input_genes = [item.upper() for item in input_genes]
        self.input_genes = input_genes

    def convert_to_Entrez(self):
        self.convert_IDs, df_convert_out = _geneplexus.initial_ID_convert(self.input_genes, self.file_loc)
        self.df_convert_out, self.table_summary, self.input_count = _geneplexus.make_validation_df(
            df_convert_out,
            self.file_loc,
        )
        return self.df_convert_out

    def set_params(self, net_type, features, GSC):
        self.net_type = net_type
        self.features = features
        self.GSC = GSC

    def get_pos_and_neg_genes(self):
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
        self.df_sim_GO, self.df_sim_Dis, self.weights_GO, self.weights_Dis = _geneplexus.make_sim_dfs(
            self.file_loc,
            self.mdl_weights,
            self.GSC,
            self.net_type,
            self.features,
        )
        return self.df_sim_GO, self.df_sim_Dis, self.weights_GO, self.weights_Dis

    def make_small_edgelist(self, num_nodes=50):
        self.df_edge, self.isolated_genes, self.df_edge_sym, self.isolated_genes_sym = _geneplexus.make_small_edgelist(
            self.file_loc,
            self.df_probs,
            self.net_type,
            num_nodes=50,
        )
        return self.df_edge, self.isolated_genes, self.df_edge_sym, self.isolated_genes_sym

    def alter_validation_df(self):
        self.df_convert_out_subset, self.positive_genes = _geneplexus.alter_validation_df(
            self.df_convert_out,
            self.table_summary,
            self.net_type,
        )
        return self.df_convert_out_subset, self.positive_genes
