import json
import os.path as osp
import shutil
import tempfile
import unittest

import numpy as np
import pytest

import geneplexus.custom


@pytest.mark.usefixtures("data")
class TestCustom(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.edgelist_loc = osp.join(pytest.HOMEDIR, "test", "custom_net.edg")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)

    @pytest.mark.order(0)
    def test_edgelist_to_nodeorder(self):
        geneplexus.custom.edgelist_to_nodeorder(
            self.edgelist_loc,
            self.tmpdir,
            "custom",
        )
        outpath = osp.join(self.tmpdir, "NodeOrder_custom.txt")
        self.assertEqual(
            np.loadtxt(outpath, dtype=str).tolist(),
            ["1213", "156", "1759", "408", "4734", "6714"],
        )

    @pytest.mark.order(1)
    def test_edgelist_to_matrix(self):
        geneplexus.custom.edgelist_to_matrix(
            self.edgelist_loc,
            self.tmpdir,
            "custom",
            "Adjacency",
        )
        outpath = osp.join(self.tmpdir, "Data_Adjacency_custom.npy")
        mat = np.load(outpath)
        self.assertEqual(
            mat.tolist(),
            [
                [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0, 0.0, 1.0],
                [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                [1.0, 1.0, 1.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            ],
        )

    @pytest.mark.order(2)
    def test_subset_GSC_to_network(self):
        geneplexus.custom.subset_GSC_to_network(
            pytest.DATADIR,
            "custom",
            "GO",
        )
        outpath = osp.join(pytest.DATADIR, "GSC_GO_custom_GoodSets.json")

        with open(outpath, "r") as f:
            goodsets = json.load(f)

        self.assertEqual(
            sorted(goodsets),
            [
                "GO:0007154",
                "GO:0007165",
                "GO:0008150",
                "GO:0009987",
                "GO:0023052",
                "GO:0050789",
                "GO:0050794",
                "GO:0050896",
                "GO:0051716",
                "GO:0065007",
            ],
        )


if __name__ == "__main__":
    unittest.main()
