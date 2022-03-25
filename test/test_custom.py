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

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)

    def test_edgelist_to_nodeorder(self):
        edgelist_loc = osp.join(pytest.HOMEDIR, "test", "custom_net.edg")
        geneplexus.custom.edgelist_to_nodeorder(
            edgelist_loc,
            self.tmpdir,
            "custom",
        )
        out_path = osp.join(self.tmpdir, "NodeOrder_custom.txt")
        self.assertEqual(
            sorted(np.loadtxt(out_path).astype(int).ravel()),
            [156, 408, 1213, 1759, 4734, 6714],
        )

    def test_edgelist_to_matrix(self):
        pass

    def test_subset_GSC_to_network(self):
        pass


if __name__ == "__main__":
    unittest.main()
