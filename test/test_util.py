import os.path as osp
import pathlib
import shutil
import tempfile
import time
import unittest

import pytest
from parameterized import parameterized

from geneplexus import config
from geneplexus import util


class TestReadGeneList(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        print(f"\nCreated temporary folder: {cls.tmpdir}")

        cls.lst = ["46", "1235", "343", "EAND"]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)
        print(f"\nDeconstructed temporary folder: {cls.tmpdir}")

    @parameterized.expand(
        [
            (",", ","),
            ("|", "|"),
            ("\t", "tab"),
            ("\n", "newline"),
        ],
    )
    def test_read(self, delimiter, delimiter_arg):
        file_path = osp.join(self.tmpdir, "lst.txt")
        with open(file_path, "w") as f:
            f.write(delimiter.join(self.lst))

        lst = util.read_gene_list(file_path, delimiter_arg)
        self.assertEqual(lst, self.lst)


class TestNetGSCGetters(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()

        fn_list = [
            "NodeOrder_customnet.txt",
            "NodeOrder_customnet2.txt",
            "GSCOriginal_customgsc.json",
            "GSCOriginal_customgsc2.json",
        ]

        for fn in fn_list:
            pathlib.Path(osp.join(cls.tmpdir, fn)).touch()

    def test_get_all_gscs(self):
        self.assertEqual(util.get_all_gscs(None), sorted(config.ALL_GSCS))
        self.assertEqual(
            util.get_all_gscs(self.tmpdir),
            sorted(config.ALL_GSCS + ["customgsc", "customgsc2"]),
        )

    def test_get_all_net_types(self):
        self.assertEqual(util.get_all_net_types(None), sorted(config.ALL_NETWORKS))
        self.assertEqual(
            util.get_all_net_types(self.tmpdir),
            sorted(config.ALL_NETWORKS + ["customnet", "customnet2"]),
        )


def test_timeout():
    @util.timeout(5)
    def wait():
        time.sleep(0.1)

    wait()

    with pytest.raises(TimeoutError):

        @util.timeout(1)
        def wait():
            time.sleep(2)

        wait()


if __name__ == "__main__":
    unittest.main()
