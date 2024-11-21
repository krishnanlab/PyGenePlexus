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
``config.yaml``               Configuration file containing the parameters used to generate the results
``df_probs.tsv``              Top predicted genes related to the input gene list.
                              (see :meth:`geneplexus.GenePlexus.fit_and_predict`)
``df_sim.tsv``                Similarity of model trained on user gene list to models trained on
                              known gene sets. (see :meth:`geneplexus.GenePlexus.make_sim_dfs`)
``df_edge.tsv``               Edgelist (Entrez ID) of subgraph induced by top predicted genes
                              (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
``df_edge_sym.tsv``           Edgelist (Symbol) of subgraph induced by top predicted genes
                              (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
``isoloated_genes.txt``       List of top predicted genes (Entrez ID) that have no edges in
                              the network. (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
``isoloated_genes_sym.txt``   List of top predicted genes (Symbol) that have no edges in
                              the network. (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
``avgps.txt``                 Cross validation evaluation of the model's ability to capture the
                              input gene list, mesured using ``log2(auprc/prior)``.
                              (see :meth:`geneplexus.GenePlexus.fit_and_predict`)
``mdl_weights.txt``           The coefficients of the trained model.
                              (see :meth:`geneplexus.GenePlexus.fit_and_predict`)
``df_convert_out.tsv``        Table showing conversion of input genes to Entrez IDs for all networks.
                              (see :meth:`geneplexus.GenePlexus.load_genes`)
``df_convert_out_subset.tsv`` Table showing conversion of input genes for network used in training.
                              (see :meth:`geneplexus.GenePlexus.alter_validation_df`)
``pos_genes_in_net.txt``      List of positive genes used in training.
                              (see :meth:`geneplexus.GenePlexus._get_pos_and_neg_genes`)
``negative_genes.txt``        List of negative genes used in training.
                              (see :meth:`geneplexus.GenePlexus._get_pos_and_neg_genes`)
``net_genes.txt``             List of all genes in the training species network.
                              (see :meth:`geneplexus.GenePlexus._get_pos_and_neg_genes`)
``neutral_gene_info.json``    Information on which genes are considered neutral (i.e. not used in training).
                              (see :meth:`geneplexus.GenePlexus._get_pos_and_neg_genes`)
``run.log``                   Run log file.
============================= ====================================================================

Full CLI options (check out with ``geneplexus --help``)

.. code-block:: text

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