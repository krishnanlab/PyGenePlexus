import os.path as osp
import shutil
import tempfile
import time
import unittest

import pytest
from parameterized import parameterized

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
