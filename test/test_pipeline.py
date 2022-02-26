import os.path as osp
import pathlib
import unittest

import pandas as pd
import pytest

from geneplexus import geneplexus


HOMEDIR = pathlib.Path(__file__).parent.parent
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
    def setUp(self):
        self.gp = geneplexus.GenePlexus(self.datadir)

    def test_filenames(self):
        for filename in FILENAMES:
            self.assertTrue(osp.isfile(osp.join(self.datadir, filename)))

    def test_pipeline(self):
        # TODO: make a helper for this
        with open(osp.join(HOMEDIR, "input_genes.txt"), "r") as f:
            input_genes = [i.strip("'") for i in f.read().split(", ")]
        self.gp.load_genes(input_genes)

        df_convert_out = self.gp.convert_to_Entrez()
        df_convert_out[["Original ID", "Entrez ID"]] = df_convert_out[["Original ID", "Entrez ID"]].astype(int)
        df_convert_out_expected = pd.read_csv(
            osp.join(HOMEDIR, "test", "expected_result", "df_convert_out.tsv"),
            sep="\t",
        )
        self.assertEqual(
            df_convert_out.values.tolist(),
            df_convert_out_expected.values.tolist(),
        )

        self.gp.set_params("BioGRID", "Embedding", "GO")
        pos_genes_in_net, negative_genes, net_genes = self.gp.get_pos_and_neg_genes()
        mdl_weights, df_probs, avgps = self.gp.fit_and_predict()


if __name__ == "__main__":
    unittest.main()
