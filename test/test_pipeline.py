import os.path as osp
import unittest

import pandas as pd
import pytest

import geneplexus


FILENAMES = [
    "CorrectionMatrix_DisGeNet_DisGeNet_BioGRID_Embedding.npy",
    "CorrectionMatrix_DisGeNet_GO_BioGRID_Embedding.npy",
    "CorrectionMatrix_GO_DisGeNet_BioGRID_Embedding.npy",
    "CorrectionMatrix_GO_GO_BioGRID_Embedding.npy",
    "CorrectionMatrixOrder_DisGeNet_BioGRID.txt",
    "CorrectionMatrixOrder_GO_BioGRID.txt",
    "Data_Embedding_BioGRID.npy",
    "Edgelist_BioGRID.edg",
    "GSC_DisGeNet_BioGRID_GoodSets.json",
    "GSC_DisGeNet_BioGRID_universe.txt",
    "GSC_GO_BioGRID_GoodSets.json",
    "GSC_GO_BioGRID_universe.txt",
    "IDconversion_Homo-sapiens_ENSG-to-Entrez.json",
    "IDconversion_Homo-sapiens_ENSP-to-Entrez.json",
    "IDconversion_Homo-sapiens_ENST-to-Entrez.json",
    "IDconversion_Homo-sapiens_Entrez-to-ENSG.json",
    "IDconversion_Homo-sapiens_Entrez-to-Name.json",
    "IDconversion_Homo-sapiens_Entrez-to-Symbol.json",
    "IDconversion_Homo-sapiens_Symbol-to-Entrez.json",
    "NodeOrder_BioGRID.txt",
    "NodeOrder_GIANT-TN.txt",
    "NodeOrder_STRING-EXP.txt",
    "NodeOrder_STRING.txt",
    "PreTrainedWeights_DisGeNet_BioGRID_Embedding.json",
    "PreTrainedWeights_GO_BioGRID_Embedding.json",
]


@pytest.mark.usefixtures("data")
def test_download_exist(caplog):
    geneplexus.download.download_select_data(
        pytest.DATADIR,
        "All",
        "BioGRID",
        "Embedding",
        ["GO", "DisGeNet"],
    )
    assert "File exists, skipping download:" in caplog.text


@pytest.mark.skip
@pytest.mark.usefixtures("data")
class TestGenePlexusPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gp = geneplexus.GenePlexus(pytest.DATADIR, "BioGRID", "Embedding", "GO")

    @pytest.mark.order(0)
    def test_filenames(self):
        for filename in FILENAMES:
            with self.subTest(filename=filename):
                self.assertTrue(osp.isfile(osp.join(pytest.DATADIR, filename)))

    @pytest.mark.order(1)
    def test_init_geneplexus(self):
        input_path = osp.join(pytest.HOMEDIR, "example", "input_genes.txt")
        input_genes = geneplexus.util.read_gene_list(input_path)
        self.gp.load_genes(input_genes)
        self.assertEqual(self.gp.input_genes, input_genes)

    @pytest.mark.order(2)
    def test_convert_to_entrez(self):
        self.gp.convert_to_Entrez()
        df_convert_out = self.gp.df_convert_out.copy()
        columns = ["Original ID", "Entrez ID"]
        df_convert_out[columns] = df_convert_out[columns].astype(int)

        path = osp.join(pytest.ANSWERDIR, "df_convert_out.tsv")
        df_convert_out_expected = pd.read_csv(path, sep="\t")
        self.assertEqual(
            df_convert_out.values.tolist(),
            df_convert_out_expected.values.tolist(),
        )

    @pytest.mark.order(3)
    def test_get_pos_and_neg_genes(self):
        self.gp.get_pos_and_neg_genes()

    @pytest.mark.order(4)
    def test_fit_and_predict(self):
        # First check if the gene IDs and the corresponding attributes are
        # aligned. Then check if the computed probabilities are close to the
        # expected results (up to third places)
        self.gp.fit_and_predict()
        df_probs = self.gp.df_probs.copy()
        df_probs["Entrez"] = df_probs["Entrez"].astype(int)

        path = osp.join(pytest.ANSWERDIR, "df_probs.tsv")
        df_probs_expected = pd.read_csv(path, sep="\t")

        # Ignore rank and prob for now as they might be susceptible to randomns
        columns = ["Entrez", "Symbol", "Name", "Known/Novel", "Class-Label"]
        self.assertEqual(
            df_probs.sort_values("Entrez")[columns].values.tolist(),
            df_probs_expected.sort_values("Entrez")[columns].values.tolist(),
        )

        # But at least the proababilities should be close
        for prob, prob_expected in zip(
            df_probs.sort_values("Entrez")["Probability"],
            df_probs_expected.sort_values("Entrez")["Probability"],
        ):
            self.assertAlmostEqual(prob, prob_expected, places=3)

    @pytest.mark.order(5)
    def test_make_sim_dfs(self):
        # First check if ID and Name are aligned. Then check if the computed
        # similarities are close to the expected results (up to third places)
        df_sim_GO, df_sim_Dis, _, _ = self.gp.make_sim_dfs()
        columns = ["ID", "Name"]

        with self.subTest("GO similarity"):
            path = osp.join(pytest.ANSWERDIR, "df_sim_GO.tsv")
            df_sim_GO_expected = pd.read_csv(path, sep="\t").fillna("NA")
            self.assertEqual(
                df_sim_GO.sort_values("ID")[columns].values.tolist(),
                df_sim_GO_expected.sort_values("ID")[columns].values.tolist(),
            )
            for sim, sim_expected in zip(
                df_sim_GO.sort_values("ID")["Similarity"],
                df_sim_GO_expected.sort_values("ID")["Similarity"],
            ):
                self.assertAlmostEqual(sim, sim_expected, places=3)

        with self.subTest("Dis similarity"):
            path = osp.join(pytest.ANSWERDIR, "df_sim_Dis.tsv")
            df_sim_Dis_expected = pd.read_csv(path, sep="\t").fillna("NA")
            self.assertEqual(
                df_sim_Dis.sort_values("ID")[columns].values.tolist(),
                df_sim_Dis_expected.sort_values("ID")[columns].values.tolist(),
            )
            for sim, sim_expected in zip(
                df_sim_Dis.sort_values("ID")["Similarity"],
                df_sim_Dis_expected.sort_values("ID")["Similarity"],
            ):
                self.assertAlmostEqual(sim, sim_expected, places=3)

    @pytest.mark.order(6)
    def test_make_small_edgelist(self):
        self.gp.make_small_edgelist(num_nodes=50)
        columns = ["Node1", "Node2"]

        with self.subTest("Edge"):
            df_edge = self.gp.df_edge.copy()
            df_edge[columns] = df_edge[columns].astype(int)
            path = osp.join(pytest.ANSWERDIR, "df_edge.tsv")
            df_edge_expected = pd.read_csv(path, sep="\t")
            self.assertEqual(
                df_edge[columns].values.tolist(),
                df_edge_expected[columns].values.tolist(),
            )
            for weight, weight_expected in zip(
                df_edge["Weight"],
                df_edge_expected["Weight"],
            ):
                self.assertAlmostEqual(weight, weight_expected)

        with self.subTest("Edge sym"):
            df_edge_sym = self.gp.df_edge_sym.copy()
            path = osp.join(pytest.ANSWERDIR, "df_edge_sym.tsv")
            df_edge_sym_expected = pd.read_csv(path, sep="\t")
            self.assertEqual(
                df_edge_sym[columns].values.tolist(),
                df_edge_sym_expected[columns].values.tolist(),
            )
            for weight, weight_expected in zip(
                df_edge_sym["Weight"],
                df_edge_sym_expected["Weight"],
            ):
                self.assertAlmostEqual(weight, weight_expected)

    @pytest.mark.order(7)
    def test_alter_validation_df(self):
        self.gp.alter_validation_df()
        df_convert_out_subset = self.gp.df_convert_out_subset.copy()
        columns = ["Original ID", "Entrez ID"]
        df_convert_out_subset[columns] = df_convert_out_subset[columns].astype(int)

        path = osp.join(pytest.ANSWERDIR, "df_convert_out_subset.tsv")
        df_convert_out_subset_expected = pd.read_csv(path, sep="\t")
        self.assertEqual(
            df_convert_out_subset.values.tolist(),
            df_convert_out_subset_expected.values.tolist(),
        )


NET_TEST_PAIRS = [
    ("BioGRID", True),
    ("STRING", True),
    ("STRING-EXP", True),
    ("GIANT-TN", True),
    ("GiANT-TN", False),
    ("TRiNG", False),
]
FEATURE_TEST_PAIRS = [
    ("Adjacency", True),
    ("Embedding", True),
    ("Influence", True),
    ("adjd", False),
    ("randomStufF", False),
]
GSC_TEST_PAIRS = [
    ("GO", True),
    ("DisGeNet", True),
    ("gO", False),
    ("CrAzyStuff", False),
]


@pytest.mark.parametrize("gsc,gsc_correct", GSC_TEST_PAIRS)
@pytest.mark.parametrize("features,features_correct", FEATURE_TEST_PAIRS)
@pytest.mark.parametrize("net_type,net_correct", NET_TEST_PAIRS)
def test_geneplexus_param(
    net_type,
    net_correct,
    features,
    features_correct,
    gsc,
    gsc_correct,
):
    if net_correct and features_correct and gsc_correct:
        geneplexus.GenePlexus(net_type=net_type, features=features, gsc=gsc)
    else:
        with pytest.raises(ValueError):
            geneplexus.GenePlexus(net_type=net_type, features=features, gsc=gsc)


if __name__ == "__main__":
    unittest.main()
