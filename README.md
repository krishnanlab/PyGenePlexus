[![Tests](https://github.com/krishnanlab/PyGenePlexus/actions/workflows/tests.yml/badge.svg)](https://github.com/krishnanlab/PyGenePlexus/actions/workflows/tests.yml)
[![Documentation Status](https://readthedocs.org/projects/pygeneplexus/badge/?version=v2.0.3)](https://pygeneplexus.readthedocs.io/en/v2.0.3/?badge=v2.0.3)
[![PyPI](https://img.shields.io/pypi/v/geneplexus)](https://pypi.org/project/geneplexus/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/geneplexus)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# PyGenePlexus [![DOI](https://zenodo.org/badge/423591778.svg)](https://zenodo.org/badge/latestdoi/423591778)

A Python package of the GenePlexus analysis pipeline.

* The [GenePlexusZoo Method Paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011773)
* The [GenePlexus Benchmarking Paper](https://academic.oup.com/bioinformatics/article/36/11/3457/5780279)
* The [Webserver](https://www.geneplexus.net/)
* The [Package Documentation](https://pygeneplexus.readthedocs.io)
* Data [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14750555.svg)](https://doi.org/10.5281/zenodo.14750555)

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
geneplexus --input_file example/input_genes.txt --output_dir example_result
```

Full CLI options (check out with ``geneplexus --help``)

```txt
Run the GenePlexus pipline on a input gene list.

options:
  -h, --help            show this help message and exit
  -i , --input_file     Input gene list (.txt) file. (default: None)
  -d , --gene_list_delimiter
                        Delimiter used in the gene list. Use 'newline' if the genes are separated
                        by new line, and use 'tab' if the genes are seperate by tabs. Other
                        generic separator are also supported, e.g. ', '. (default: newline)
  -dd , --data_dir      Directory in which the data are stored, if set to None, then use the
                        default data directory ~/.data/geneplexus (default: None)
  -n , --network        Network to use. The choices are: {BioGRID, STRING, IMP} (default: STRING)
  -f , --feature        Types of feature to use. The choices are: {SixSpeciesN2V} (default:
                        SixSpeciesN2V)
  -s1 , --sp_trn        Species of training data The choices are: {Human, Mouse, Fly, Worm,
                        Zebrafish, Yeast} (default: Human)
  -s2 , --sp_res        Species of results data The choices are: {Human, Mouse, Fly, Worm,
                        Zebrafish, Yeast} (default: Mouse)
  -g1 , --gsc_trn       Geneset collection used to generate negatives. The choices are: {GO,
                        Monarch, Mondo, Combined} (default: GO)
  -g2 , --gsc_res       Geneset collection used for model similarities. The choices are: {GO,
                        Monarch, Mondo, Combined} (default: GO)
  -s , --small_edgelist_num_nodes
                        Number of nodes in the small edgelist. (default: 50)
  -od , --output_dir    Output directory with respect to the repo root directory. (default:
                        result/)
  -l , --log_level      Logging level. The choices are: {CRITICAL, ERROR, WARNING, INFO, DEBUG}
                        (default: INFO)
  -ad, --auto_download_off
                        Turns off autodownloader which is on by default. (default: False)
  -q, --quiet           Suppress log messages (same as setting log_level to CRITICAL). (default:
                        False)
  -z, --zip-output      If set, then compress the output directory into a Zip file. (default:
                        False)
  --clear-data          Clear data directory and exit. (default: False)
  --overwrite           Overwrite existing result directory if set. (default: False)
  --skip-mdl-sim        Skip model similarity computation (default: False)
  --skip-sm-edgelist    Skip making small edgelist. (default: False)
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

1. Install doc dependencies

```bash
pip install -r docs/requirements.txt
```

2. Build

```bash
cd docs
make html
```

3. Open doc

```bash
open build/html/index.html
```
