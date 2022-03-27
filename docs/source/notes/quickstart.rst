Quick start
===========

PyGenePlexus CLI
----------------

PyGenePlexus provides a command line interface to run the full [GenePlexus]_
pipeline on a user defined geneset (a text file with gene ids seprated by line).

.. code-block:: bash

   geneplexus  --input_file my_gene_list.txt --output_dir result --data_dir data

The command above read the gene list file ``my_gene_list.txt``, download necessary
data files and save under ``data``, and generates the following files under ``result/``

============================= ====================================================================
``df_convert_out.tsv``        Entrez gene ID conversion result table
                              (see :meth:`geneplexus.GenePlexus.convert_to_Entrez`)
``df_edge.tsv``               Edgelist (gene Entrez ID) of subgraph induced by top predicted genes
                              (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
``df_edge_sym.tsv``           Edgelist (gene symbol) of subgraph induced by top predicted genes
                              (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
``df_probs.tsv``              Top predicted genes related to the input gene list
                              (see :meth:`geneplexus.GenePlexus.fit_and_predict`)
``df_sim_GO.tsv``             Model similarity with those from the GO terms
                              (see :meth:`geneplexus.GenePlexus.make_sim_dfs`)
``df_sim_Dis.tsv``            Model similarity with those from the DisGeNet terms
                              (see :meth:`geneplexus.GenePlexus.make_sim_dfs`)
``df_convert_out_subset.tsv`` See :meth:`geneplexus.GenePlexus.alter_validation_df`
============================= ====================================================================

Full CLI options (check out with ``geneplexus --help``)

.. code-block:: text

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
                            also supported, e.g. ', '. (default: newline)
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


PyGenePlexus API
----------------

Downloading datasets
^^^^^^^^^^^^^^^^^^^^

Running the PyGenePlexus pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
