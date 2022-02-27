import os.path as osp
import pathlib
import unittest

import pandas as pd
import pytest

from geneplexus import geneplexus


HOMEDIR = pathlib.Path(__file__).parent.parent
ANSWERDIR = osp.join(HOMEDIR, "test", "expected_result")
FILENAMES = [
    "CorrectionMatrix_DisGeNet_DisGeNet_BioGRID_Embedding.npy",
    "CorrectionMatrix_DisGeNet_GO_BioGRID_Embedding.npy",
    "CorrectionMatrix_GO_DisGeNet_BioGRID_Embedding.npy",
    "CorrectionMatrix_GO_GO_BioGRID_Embedding.npy",
    "CorrectionMatrixOrder_DisGeNet_BioGRID.txt",
    "CorrectionMatrixOrder_GO_BioGRID.txt",
    "Data_Embedding_BioGRID.npy",
    "Edgelist_BioGRID.edg",
    "GSC_DisGeNet_BioGRID_GoodSets.pickle",
    "GSC_DisGeNet_BioGRID_universe.txt",
    "GSC_GO_BioGRID_GoodSets.pickle",
    "GSC_GO_BioGRID_universe.txt",
    "IDconversion_Homo-sapiens_ENSG-to-Entrez.pickle",
    "IDconversion_Homo-sapiens_ENSP-to-Entrez.pickle",
    "IDconversion_Homo-sapiens_ENST-to-Entrez.pickle",
    "IDconversion_Homo-sapiens_Entrez-to-ENSG.pickle",
    "IDconversion_Homo-sapiens_Entrez-to-Name.pickle",
    "IDconversion_Homo-sapiens_Entrez-to-Symbol.pickle",
    "IDconversion_Homo-sapiens_Symbol-to-Entrez.pickle",
    "NodeOrder_BioGRID.txt",
    "NodeOrder_GIANT-TN.txt",
    "NodeOrder_STRING-EXP.txt",
    "NodeOrder_STRING.txt",
    "PreTrainedWeights_DisGeNet_BioGRID_Embedding.pickle",
    "PreTrainedWeights_GO_BioGRID_Embedding.pickle",
]


@pytest.fixture(scope="class")
def data(request):
    datadir = request.config.cache.makedir("download")
    geneplexus.download_select_data(
        datadir,
        "All",
        "BioGRID",
        "Embedding",
        ["GO", "DisGeNet"],
    )
    request.cls.datadir = datadir


@pytest.mark.usefixtures("data")
class TestGenePlexusPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gp = geneplexus.GenePlexus(HOMEDIR)
        cls.gp.set_params("BioGRID", "Embedding", "GO")

    @pytest.mark.order(0)
    def test_filenames(self):
        for filename in FILENAMES:
            self.assertTrue(osp.isfile(osp.join(self.datadir, filename)))

    @pytest.mark.order(1)
    def test_init_geneplexus(self):
        self.gp.file_loc = self.datadir
        input_genes_path = osp.join(HOMEDIR, "input_genes.txt")
        input_genes = geneplexus.util.read_gene_list(input_genes_path)
        self.gp.load_genes(input_genes)
        self.assertEqual(self.gp.input_genes, input_genes)

    @pytest.mark.order(2)
    def test_convert_to_entrez(self):
        self.gp.convert_to_Entrez()
        df_convert_out = self.gp.df_convert_out.copy()
        columns = ["Original ID", "Entrez ID"]
        df_convert_out[columns] = df_convert_out[columns].astype(int)

        path = osp.join(ANSWERDIR, "df_convert_out.tsv")
        df_convert_out_expected = pd.read_csv(path, sep="\t")
        self.assertEqual(
            df_convert_out.values.tolist(),
            df_convert_out_expected.values.tolist(),
        )

    @pytest.mark.order(3)
    def test_get_pos_and_neg_genes(self):
        self.gp.get_pos_and_neg_genes()

    @pytest.mark.xfail(reason="Randomness, need to set random state.")
    @pytest.mark.order(4)
    def test_fit_and_predict(self):
        self.gp.fit_and_predict()
        df_probs = self.gp.df_probs.copy()
        df_probs["Entrez"] = df_probs["Entrez"].astype(int)

        path = osp.join(ANSWERDIR, "df_probs.tsv")
        df_probs_expected = pd.read_csv(path, sep="\t")
        columns = ["Entrez", "Symbol", "Name", "Known/Novel", "Class-Label", "Rank"]
        self.assertEqual(
            df_probs[columns].values.tolist(),
            df_probs_expected[columns].values.tolist(),
        )

        for prob, prob_expected in zip(
            df_probs["Probability"],
            df_probs_expected["Probability"],
        ):
            self.assertAlmostEqual(prob, prob_expected)

    @pytest.mark.xfail(reason="Randomness, need to set random state.")
    @pytest.mark.order(5)
    def test_make_sim_dfs(self):
        df_sim_GO, df_sim_Dis, _, _ = self.gp.make_sim_dfs()
        columns = ["ID", "Rank"]

        with self.subTest("GO similarity"):
            path = osp.join(ANSWERDIR, "df_sim_GO.tsv")
            df_sim_GO_expected = pd.read_csv(path, sep="\t")
            self.assertEqual(
                df_sim_GO[columns].values.tolist(),
                df_sim_GO_expected[columns].values.tolist(),
            )
            for sim, sim_expected in zip(
                df_sim_GO["Similarity"],
                df_sim_GO_expected["Similarity"],
            ):
                self.assertAlmostEqual(sim, sim_expected)

        with self.subTest("Dis similarity"):
            path = osp.join(ANSWERDIR, "df_sim_Dis.tsv")
            df_sim_Dis_expected = pd.read_csv(path, sep="\t")
            self.assertEqual(
                df_sim_Dis[columns].values.tolist(),
                df_sim_Dis_expected[columns].values.tolist(),
            )
            for sim, sim_expected in zip(
                df_sim_Dis["Similarity"],
                df_sim_Dis_expected["Similarity"],
            ):
                self.assertAlmostEqual(sim, sim_expected)

    @pytest.mark.order(6)
    def test_make_small_edgelist(self):
        self.gp.make_small_edgelist(num_nodes=50)
        columns = ["Node1", "Node2"]

        with self.subTest("Edge"):
            df_edge = self.gp.df_edge.copy()
            df_edge[columns] = df_edge[columns].astype(int)
            path = osp.join(ANSWERDIR, "df_edge.tsv")
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
            path = osp.join(ANSWERDIR, "df_edge_sym.tsv")
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

        path = osp.join(ANSWERDIR, "df_convert_out_subset.tsv")
        df_convert_out_subset_expected = pd.read_csv(path, sep="\t")
        self.assertEqual(
            df_convert_out_subset.values.tolist(),
            df_convert_out_subset_expected.values.tolist(),
        )


if __name__ == "__main__":
    unittest.main()
