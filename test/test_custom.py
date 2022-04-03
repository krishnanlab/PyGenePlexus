import json
import os
import os.path as osp
import pathlib
import shutil
import tempfile
import unittest

import numpy as np
import pytest
from parameterized import parameterized

import geneplexus


TESTDIR = osp.join(pathlib.Path(__file__).absolute().parent)
EDGELIST_UNWEIGHTED_LOC = osp.join(TESTDIR, "custom_net.edg")
EDGELIST_WEIGHTED_LOC = osp.join(TESTDIR, "custom_net_weighted.edg")

NODEORDER = ["1213", "156", "1759", "408", "4734", "6714"]
ADJMAT_UNWEIGHTED = [
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 1.0, 0.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
]
ADJMAT_WEIGHTED = [
    [0.0, 0.0, 0.0, 0.3, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.4, 0.0, 0.6],
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    [0.3, 0.4, 1.0, 0.0, 0.6, 0.0],
    [0.0, 0.0, 0.0, 0.6, 0.0, 0.0],
    [0.0, 0.6, 0.0, 0.0, 0.0, 0.0],
]


@pytest.mark.usefixtures("data")
class TestCustom(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.nodeorder_path = osp.join(pytest.DATADIR, "NodeOrder_custom.txt")
        np.savetxt(cls.nodeorder_path, NODEORDER, fmt="%s")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)
        os.remove(cls.nodeorder_path)

    @parameterized.expand(
        [
            (EDGELIST_UNWEIGHTED_LOC,),
            (EDGELIST_WEIGHTED_LOC,),
        ],
    )
    def test_edgelist_to_nodeorder(self, edgelist_loc):
        geneplexus.custom.edgelist_to_nodeorder(edgelist_loc, self.tmpdir, "custom")
        outpath = osp.join(pytest.DATADIR, "NodeOrder_custom.txt")
        self.assertEqual(np.loadtxt(outpath, dtype=str).tolist(), NODEORDER)

    @parameterized.expand(
        [
            (EDGELIST_UNWEIGHTED_LOC, ADJMAT_UNWEIGHTED),
            (EDGELIST_WEIGHTED_LOC, ADJMAT_WEIGHTED),
        ],
    )
    def test_edgelist_to_matrix(self, edgelist_loc, adjmat):
        geneplexus.custom.edgelist_to_matrix(
            edgelist_loc,
            pytest.DATADIR,
            "custom",
            "Adjacency",
        )
        outpath = osp.join(pytest.DATADIR, "Data_Adjacency_custom.npy")
        self.assertEqual(np.load(outpath).tolist(), adjmat)

        # TODO: Move this to a test func
        geneplexus.GenePlexus(pytest.DATADIR, "custom", "Adjacency")

    def test_subset_gsc_to_network(self):
        geneplexus.custom.subset_gsc_to_network(
            pytest.DATADIR,
            "custom",
            "GO",
            min_size=6,
        )
        outpath = osp.join(pytest.DATADIR, "GSC_GO_custom_GoodSets.json")

        with open(outpath, "r") as f:
            goodsets = json.load(f)

        self.assertEqual(
            sorted(goodsets),
            [
                "GO:0006810",
                "GO:0007154",
                "GO:0007165",
                "GO:0008150",
                "GO:0009987",
                "GO:0016192",
                "GO:0023052",
                "GO:0050789",
                "GO:0050794",
                "GO:0050896",
                "GO:0051179",
                "GO:0051234",
                "GO:0051716",
                "GO:0065007",
            ],
        )


if __name__ == "__main__":
    unittest.main()
