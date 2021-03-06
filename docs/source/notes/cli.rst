PyGenePlexus CLI
================

PyGenePlexus provides a command line interface to run the full [GenePlexus]_
pipeline on a user defined geneset (a text file with gene IDs seprated by line).

.. code-block:: bash

   geneplexus  --input_file my_gene_list.txt --output_dir my_result --data_dir my_data

The command above reads the gene list file ``my_gene_list.txt``, downloads the necessary
data files and saves them to the directory ``my_data/``. If ``--data_dir`` is not supplied,
the data files will be saved under ``~/.data/geneplexus/`` by default. Finally, all
output files will be saved under ``my_result/``.

.. note::

    If the direcory ``my_result/`` already exists, the program will try to append
    a number, e.g., ``my_result_1/``, to prevent overwriting. If you would like
    to overwrite, you can do so by specifying the ``--overwrite`` CLI option.

The output files contain the following files.

============================= ====================================================================
``config.yaml``               Configuration file containing the parameters used to generate the results and the input
                              gene list.
``cross_validation.txt``      Cross validation evaluation of the model's ability to capture the
                              input gene list, mesured using ``log2(auprc/prior)``
``df_convert_out.tsv``        Table showing conversion of input genes to Entrez IDs
                              (see :meth:`geneplexus.GenePlexus.load_genes`)
``df_edge.tsv``               Edgelist (gene Entrez ID) of subgraph induced by top predicted genes
                              (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
``df_edge_sym.tsv``           Edgelist (gene symbol) of subgraph induced by top predicted genes
                              (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
``df_probs.tsv``              Top predicted genes related to the input gene list
                              (see :meth:`geneplexus.GenePlexus.fit_and_predict`)
``df_sim_GO.tsv``             Similarity of model trained on user gene list to models trained on
                              known GO gene sets. (see :meth:`geneplexus.GenePlexus.make_sim_dfs`)
``df_sim_Dis.tsv``            Similarity of model trained on user gene list to models trained on
                              known DisGeNet gene sets (see :meth:`geneplexus.GenePlexus.make_sim_dfs`)
``df_convert_out_subset.tsv`` See :meth:`geneplexus.GenePlexus.alter_validation_df`
``run.log``                   Run log file.
============================= ====================================================================

Full CLI options (check out with ``geneplexus --help``)

.. code-block:: text

    Run the GenePlexus pipline on a input gene list.

    optional arguments:
      -h, --help            show this help message and exit
      -i , --input_file     Input gene list (.txt) file (one gene per line). (default: None)
      -d , --gene_list_delimiter
                            Delimiter used in the gene list. Use 'newline' if the genes are separated
                            by new line, and use 'tab' if the genes are seperate by tabs. Other
                            generic separator are also supported, e.g. ', '. (default: newline)
      -n , --network        Network to use. {format_choices(config.ALL_NETWORKS)} (default: STRING)
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
      -q, --quiet           Suppress log messages (same as setting log_level to CRITICAL). (default:
                            False)
      -z, --zip-output      If set, then compress the output directory into a Zip file. (default:
                            False)
      --clear-data          Clear data directory and exit. (default: False)
      --overwrite           Overwrite existing result directory if set. (default: False)
      --skip-mdl-sim        Skip model similarity computation. This computation is not yet available
                            when using custom networks due to the lack of pretrained models for
                            comparison. (default: False)
