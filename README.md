# GeneplexusPublic
The point of the repo is to a get code base that can be imported by ADS folks without them having to alter any of the python code

To see how this can be run see the `example_run.py` script

# Quick start

## Installation

Install the ``GenePlexus`` analysis package via ``pip``.

```bash
pip install .

# Alternatively, install in editable mode
pip install -e .
```

## Run GenePlexus pipline

```bash
geneplexus --input_file input_genes.txt --output_dir result/ --daat_dir data/
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
```
