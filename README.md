[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6369426.svg)](https://doi.org/10.5281/zenodo.6369426)
[![Documentation Status](https://readthedocs.org/projects/pygeneplexus/badge/?version=main)](https://pygeneplexus.readthedocs.io/en/main/?badge=main)
[![Tests](https://github.com/krishnanlab/PyGenePlexus/actions/workflows/tests.yml/badge.svg)](https://github.com/krishnanlab/PyGenePlexus/actions/workflows/tests.yml)

# GeneplexusPublic

A Python package of the GenePlexus analysis pipeline.

* The [GenePlexus paper](https://academic.oup.com/bioinformatics/article/36/11/3457/5780279)
* The [repository](https://github.com/krishnanlab/GenePlexus) for reproducing the experiments
* The [webserver](https://www.geneplexus.net/)

# Quick start

## Installation

Install the ``GenePlexus`` package via ``pip``.

```bash
pip install .
```

## Run GenePlexus pipline

### Example script

See `example/example_run.py` for example usage of the API.

### Command-line interface

```bash
geneplexus --input_file example/input_genes.txt --output_dir result/ --daat_dir data/
```

Full options (check out with ``geneplexus --help``)

```txt
Run the GenePlexus pipline on a input gene list.

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE
                        Input gene list (.txt) file (one gene per line).
                        (default: None)
  --gene_list_delimiter GENE_LIST_DELIMITER
                        Delimiter used in the gene list. Use 'newline' if the
                        genes are separated by new line, and use 'tab' if the
                        genes are seperate by tabs. Other generic separator are
                        also supported, e.g. ', '. (default: , )
  --network {BioGRID,STRING,STRING-EXP,GIANT-TN}
                        Network to use for generating features.
                        (default: BioGRID)
  --feature {Adjacency,Embedding,Influence}
                        Types of feature to use. (default: Embedding)
  --GSC {GO,DisGeNet}   Geneset collection used to generate negatives and the
                        model similarities. (default: GO)
  --small_edgelist_num_nodes SMALL_EDGELIST_NUM_NODES
                        Number of nodes in the small edgelist. (default: 50)
  --data_dir DATA_DIR   Directory in which the data are stored.
                        (default: data/)
  --output_dir OUTPUT_DIR
                        Output directory with respect to the repo root
                        directory. (default: result/)
  --zip_output          If set, then compress the output directory into a Tar
                        Gz file. (default: False)
```

# Dev

## Installation

Install the PyGenePlexus package in editable mode with dev dependencies

```bash
pip install -e ."[dev]"
```

## Testing

Run the default test suite

```bash
pytest test/
```

By default, test data will be cached. Thus, after the first test run, data redownload will not be tested. To force redownload, specify the ``--cache-clear`` option

```bash
pytest test/ --cache-clear
```

## Building Documentation

1. Install doc dependencies ``pip install -r docs/requirements.txt``
2. Build
```bash
cd docs
make html
```
3. Open doc ``open build/html/index.html``
