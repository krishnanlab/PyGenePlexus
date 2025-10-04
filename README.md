[![Tests](https://github.com/krishnanlab/PyGenePlexus/actions/workflows/tests.yml/badge.svg)](https://github.com/krishnanlab/PyGenePlexus/actions/workflows/tests.yml)
[![Documentation Status](https://readthedocs.org/projects/pygeneplexus/badge/?version=v2.0.4)](https://pygeneplexus.readthedocs.io/en/v2.0.4/?badge=v2.0.4)
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
