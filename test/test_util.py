import os.path as osp
import shutil
import tempfile
import unittest

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

    def test_read_csv(self):
        file_path = osp.join(self.tmpdir, "lst.txt")
        with open(file_path, "w") as f:
            f.write(",".join(self.lst))

        lst = util.read_gene_list(file_path, ",")
        self.assertEqual(lst, self.lst)

    def test_read_tsv(self):
        file_path = osp.join(self.tmpdir, "lst.txt")
        with open(file_path, "w") as f:
            f.write("\t".join(self.lst))

        lst = util.read_gene_list(file_path, "tab")
        self.assertEqual(lst, self.lst)

    def test_read_linesep(self):
        file_path = osp.join(self.tmpdir, "lst.txt")
        with open(file_path, "w") as f:
            f.write("\n".join(self.lst))

        lst = util.read_gene_list(file_path, "newline")
        self.assertEqual(lst, self.lst)

    def test_read_pipesep(self):
        file_path = osp.join(self.tmpdir, "lst.txt")
        with open(file_path, "w") as f:
            f.write("|".join(self.lst))

        lst = util.read_gene_list(file_path, "|")
        self.assertEqual(lst, self.lst)


if __name__ == "__main__":
    unittest.main()
