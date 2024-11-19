import os.path as osp
import shutil
import tempfile
import unittest

import pandas as pd
import pytest
import yaml

import geneplexus


@pytest.mark.usefixtures("data")
class TestGenePlexusPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gp = geneplexus.GenePlexus(
            file_loc=pytest.DATADIR,
            net_type="STRING",
            features="SixSpeciesN2V",
            sp_trn="Human",
            sp_res="Mouse",
            gsc_trn="Combined",
            gsc_res="Combined",
            auto_download=False,
            log_level="INFO",
        )
        cls.tmpdir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)

    @pytest.mark.order(0)
    def test_filenames(self):
        for filename in pytest.FILENAMES:
            with self.subTest(filename=filename):
                self.assertTrue(osp.isfile(osp.join(pytest.DATADIR, filename)))

    @pytest.mark.order(1)
    def test_load_genes(self):
        input_genes = geneplexus.util.read_gene_list(pytest.GENELIST_PATH)
        self.gp.load_genes(input_genes)
        self.assertEqual(self.gp.input_genes, input_genes)

    @pytest.mark.order(2)
    def test_dump_config(self):
        self.gp.dump_config(self.tmpdir)
        with open(osp.join(pytest.ANSWERDIR, "config.yaml")) as f1, open(
            osp.join(self.tmpdir, "config.yaml"),
        ) as f2:
            cfg1, cfg2 = yaml.load(f1, yaml.Loader), yaml.load(f2, yaml.Loader)
        for param in [
            "net_type",
            "features",
            "sp_trn",
            "sp_res",
            "gsc_trn",
            "gsc_res",
            "auto_download",
            "log_level",
            "input_genes",
        ]:
            self.assertEqual(cfg1[param], cfg2[param])

    @pytest.mark.order(2)
    def test_df_convert_out(self):
        df_convert_out = self.gp.df_convert_out.copy()
        columns = ["Original ID", "Entrez ID"]
        df_convert_out[columns] = df_convert_out[columns].astype(str)

        path = osp.join(pytest.ANSWERDIR, "df_convert_out.tsv")
        df_convert_out_expected = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(
            df_convert_out[df_convert_out_expected.columns].values.tolist(),
            df_convert_out_expected.values.tolist(),
        )

    @pytest.mark.order(3)
    def test_fit_and_predict(self):
        # First check if the gene IDs and the corresponding attributes are
        # aligned. Then check if the computed probabilities are close to the
        # expected results (up to third places)
        self.gp.fit_and_predict()
        df_probs = self.gp.df_probs.copy()
        df_probs["Entrez"] = df_probs["Entrez"].astype(str)

        path = osp.join(pytest.ANSWERDIR, "df_probs.tsv")
        df_probs_expected = pd.read_csv(path, sep="\t", keep_default_na=False)
        df_probs_expected["Entrez"] = df_probs_expected["Entrez"].astype(str)

        # Ignoring columns that might be susceptible to randomns
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

    @pytest.mark.order(4)
    def test_make_sim_df(self):
        # First check if Task and ID and Name are aligned. Then check if the computed
        # similarities are close to the expected results (up to third places)
        df_sim_GO, weights_dict = self.gp.make_sim_dfs()
        columns = ["Task", "ID", "Name"]

        path = osp.join(pytest.ANSWERDIR, "df_sim.tsv")
        df_sim_GO_expected = pd.read_csv(path, sep="\t", keep_default_na=False)
        self.assertEqual(
            df_sim_GO.sort_values("ID")[columns].values.tolist(),
            df_sim_GO_expected.sort_values("ID")[columns].values.tolist(),
        )

        for sim, sim_expected in zip(
            df_sim_GO.sort_values("ID")["Similarity"],
            df_sim_GO_expected.sort_values("ID")["Similarity"],
        ):
            self.assertAlmostEqual(sim, sim_expected, places=3)

    @pytest.mark.order(5)
    def test_make_small_edgelist(self):
        self.gp.make_small_edgelist(num_nodes=50)
        columns = ["Node1", "Node2"]

        with self.subTest("Edge"):
            df_edge = self.gp.df_edge.copy()
            df_edge[columns] = df_edge[columns].astype(str)
            path = osp.join(pytest.ANSWERDIR, "df_edge.tsv")
            df_edge_expected = pd.read_csv(path, sep="\t", keep_default_na=False)
            df_edge_expected[columns] = df_edge_expected[columns].astype(str)
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
            df_edge_sym_expected = pd.read_csv(path, sep="\t", keep_default_na=False)
            self.assertEqual(
                df_edge_sym[columns].values.tolist(),
                df_edge_sym_expected[columns].values.tolist(),
            )
            for weight, weight_expected in zip(
                df_edge_sym["Weight"],
                df_edge_sym_expected["Weight"],
            ):
                self.assertAlmostEqual(weight, weight_expected)

    @pytest.mark.order(6)
    def test_alter_validation_df(self):
        self.gp.alter_validation_df()
        df_convert_out_subset = self.gp.df_convert_out_subset.copy()
        columns = ["Original ID", "Entrez ID"]
        df_convert_out_subset[columns] = df_convert_out_subset[columns].astype(str)

        path = osp.join(pytest.ANSWERDIR, "df_convert_out_subset.tsv")
        df_convert_out_subset_expected = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(
            df_convert_out_subset.values.tolist(),
            df_convert_out_subset_expected.values.tolist(),
        )


if __name__ == "__main__":
    unittest.main()
