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
from geneplexus.exception import CustomDataError


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
        cls.adj_path = osp.join(pytest.DATADIR, "Data_Adjacency_custom.npy")
        cls.gsc_path = osp.join(pytest.DATADIR, "GSC_GO_custom_GoodSets.json")
        np.savetxt(cls.nodeorder_path, NODEORDER, fmt="%s")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)
        os.remove(cls.nodeorder_path)
        os.remove(cls.adj_path)
        os.remove(cls.gsc_path)

    @parameterized.expand(
        [
            (EDGELIST_UNWEIGHTED_LOC,),
            (EDGELIST_WEIGHTED_LOC,),
        ],
    )
    def test_edgelist_to_nodeorder(self, edgelist_loc):
        geneplexus.custom.edgelist_to_nodeorder(edgelist_loc, self.tmpdir, "custom")
        self.assertEqual(np.loadtxt(self.nodeorder_path, dtype=str).tolist(), NODEORDER)

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
        self.assertEqual(np.load(self.adj_path).tolist(), adjmat)

    def test_subset_gsc_to_network(self):
        geneplexus.custom.subset_gsc_to_network(
            pytest.DATADIR,
            "custom",
            "GO",
            min_size=6,
        )

        with open(self.gsc_path, "r") as f:
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


class TestCustomGenePlexus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()

        fn_list = [
            "NodeOrder_customnet.txt",
            "Data_Adjacency_customnet.npy",
            "GSC_GO_customnet_GoodSets.json",
            "GSC_GO_customnet_universe.txt",
            "GSCOriginal_customgsc.json",
            "GSC_customgsc_customnet_GoodSets.json",
            "GSC_customgsc_customnet_universe.txt",
        ]

        for fn in fn_list:
            pathlib.Path(osp.join(cls.tmpdir, fn)).touch()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)

    def test_custom_geneplexus_init_customnet(self):
        geneplexus.GenePlexus(self.tmpdir, "customnet", "Adjacency", "GO")

    def test_custom_geneplexus_init_customgsc(self):
        geneplexus.GenePlexus(self.tmpdir, "customnet", "Adjacency", "customgsc")

    def test_custom_geneplexus_init_fail_net(self):
        with self.assertRaises(ValueError) as e:
            geneplexus.GenePlexus(self.tmpdir, net_type="customnet2")
        self.assertEqual(
            str(e.exception),
            "Unexpected network 'customnet2', available choices are "
            "['BioGRID', 'GIANT-TN', 'STRING', 'STRING-EXP', 'customnet']",
        )

    def test_custom_geneplexus_init_fail_gsc(self):
        with self.assertRaises(ValueError) as e:
            geneplexus.GenePlexus(self.tmpdir, gsc="customgsc2")
        self.assertEqual(
            str(e.exception),
            "Unexpected GSC 'customgsc2', available choices are ['DisGeNet', 'GO', 'customgsc']",
        )

    @parameterized.expand(
        [
            ("Influence", "GO"),  # Influence feature not set
            ("Influence", "DisGeNet"),  # unknown gsc DisGeNet and above
            ("Adjacency", "DisGeNet"),  # unknown gsc DisGeNet
        ],
    )
    def test_custom_geneplexus_init_fail_custom(self, features, gsc):
        with self.assertRaises(CustomDataError):
            geneplexus.GenePlexus(self.tmpdir, "customnet", features, gsc)


if __name__ == "__main__":
    unittest.main()
