[![Tests](https://github.com/krishnanlab/PyGenePlexus/actions/workflows/tests.yml/badge.svg)](https://github.com/krishnanlab/PyGenePlexus/actions/workflows/tests.yml)
[![Documentation Status](https://readthedocs.org/projects/pygeneplexus/badge/?version=main)](https://pygeneplexus.readthedocs.io/en/main/?badge=main)
[![PyPI](https://img.shields.io/pypi/v/geneplexus)](https://pypi.org/project/geneplexus/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/geneplexus)
[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# PyGeneplexus

A Python package of the GenePlexus analysis pipeline.

* The [GenePlexus paper](https://academic.oup.com/bioinformatics/article/36/11/3457/5780279)
* The [repository](https://github.com/krishnanlab/GenePlexus) for reproducing the experiments
* The [webserver](https://www.geneplexus.net/)
* Data [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6383205.svg)](https://doi.org/10.5281/zenodo.6383205)

# Quick start

## Installation

Install the ``GenePlexus`` package via ``pip``.

```bash
pip install geneplexus
```

## Run GenePlexus pipline

### Example script

See `example/example_run.py` for example usage of the API.

### Command-line interface

```bash
geneplexus --input_file example/input_genes.txt --output_dir result/ --data_dir data/
```

Full CLI options (check out with ``geneplexus --help``)

```txt
Run the GenePlexus pipline on a input gene list.

optional arguments:
  -h, --help            show this help message and exit
  -i , --input_file     Input gene list (.txt) file (one gene per line). (default: None)
  -d , --gene_list_delimiter
                        Delimiter used in the gene list. Use 'newline' if the genes are separated
                        by new line, and use 'tab' if the genes are seperate by tabs. Other
                        generic separator are also supported, e.g. ', '. (default: newline)
  -n , --network        Network to use. {format_choices(config.ALL_NETWORKS)} (default: BioGRID)
  -f , --feature        Types of feature to use. The choices are: {Adjacency, Embedding,
                        Influence} (default: Embedding)
  -g , --gsc            Geneset collection used to generate negatives and the modelsimilarities.
                        The choices are: {GO, DisGeNet} (default: GO)
  -s , --small_edgelist_num_nodes
                        Number of nodes in the small edgelist. (default: 50)
  -dd , --data_dir      Directory in which the data are stored, if set to None, then use the
                        default data directory ~/.data/geneplexus (default: None)
  -od , --output_dir    Output directory with respect to the repo root directory. (default:
                        result/)
  -l , --log_level      Logging level. The choices are: {CRITICAL, ERROR, WARNING, INFO, DEBUG}
                        (default: INFO)
  -q, --quiet           Suppress log messages (same as setting lov_level to CRITICAL). (default:
                        False)
  -z, --zip-output      If set, then compress the output directory into a Zip file. (default:
                        False)
  --clear-data          Clear data directory and exit. (default: False)
  --overwrite           Overwrite existing result directory if set. (default: False)
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
