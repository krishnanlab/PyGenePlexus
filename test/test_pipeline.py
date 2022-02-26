import os.path as osp
import pathlib
import shutil
import tempfile
import unittest

import pandas as pd
import pytest

from geneplexus import geneplexus


HOMEDIR = pathlib.Path(__file__).parent.parent
FILENAMES = [
    "CorrectionMatrix_GO_GO_BioGRID_Embedding.npy",
    "CorrectionMatrixOrder_GO_BioGRID.txt",
    "Data_Embedding_BioGRID.npy",
    "Edgelist_BioGRID.edg",
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
    "PreTrainedWeights_GO_BioGRID_Embedding.pickle",
]


class TestGenePlexus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tempdir)

    @pytest.mark.dependency()
    def test_download(self):
        geneplexus.download_select_data(
            self.tempdir,
            tasks="All",
            networks="BioGRID",
            features="Embedding",
            GSCs="GO",
        )

        for filename in FILENAMES:
            self.assertTrue(osp.isfile(osp.join(self.tempdir, filename)))

    @pytest.mark.dependency(depends=["TestGenePlexus::test_download"])
    def test_pipeline(self):
        myclass = geneplexus.GenePlexus(self.tempdir)

        # TODO: make a helper for this
        with open(osp.join(HOMEDIR, "input_genes.txt"), "r") as f:
            input_genes = [i.strip("'") for i in f.read().split(", ")]
        myclass.load_genes(input_genes)

        df_convert_out = myclass.convert_to_Entrez()
        df_convert_out[["Original ID", "Entrez ID"]] = df_convert_out[["Original ID", "Entrez ID"]].astype(int)
        df_convert_out_expected = pd.read_csv(
            osp.join(HOMEDIR, "test", "expected_result", "df_convert_out.tsv"),
            sep="\t",
        )
        self.assertEqual(
            df_convert_out.values.tolist(),
            df_convert_out_expected.values.tolist(),
        )

        myclass.set_params("BioGRID", "Embedding", "GO")
        pos_genes_in_net, negative_genes, net_genes = myclass.get_pos_and_neg_genes()
        mdl_weights, df_probs, avgps = myclass.fit_and_predict()


if __name__ == "__main__":
    unittest.main()
