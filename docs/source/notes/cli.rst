.. _cli:

PyGenePlexus CLI
================

PyGenePlexus provides a command line interface to run the full [GenePlexus]_
pipeline on a user defined geneset (a text file with gene IDs seprated by line).

.. code-block:: bash

   geneplexus  --input_file my_gene_list.txt --output_dir my_result --file_loc my_data

The command above reads the gene list file ``my_gene_list.txt``, downloads the necessary
data files and saves them to the directory ``my_data/``. If ``--file_loc`` is not supplied,
the data files will be saved under ``~/.data/geneplexus/`` by default. Finally, all
output files will be saved under ``my_result/``.

.. note::

    If the direcory ``my_result/`` already exists, the program will try to append
    a number, e.g., ``my_result_1/``, to prevent overwriting. If you would like
    to overwrite, you can do so by specifying the ``--overwrite`` CLI option.


Full CLI options (check out with ``geneplexus --help``)

.. code-block:: text

	Run the GenePlexus pipline on a input gene list.

	options:
	  -h, --help            show this help message and exit
	  -i , --input_file     Input gene list file (eg. (.txt file)). (default: None)
	  -d , --gene_list_delimiter
	                        Delimiter used in the gene list. Use 'newline' if the genes are separated
	                        by new line, and use 'tab' if the genes are seperate by tabs. If not
	                        newline or tab, will use argument directly, so /t, /n, , (default:
	                        newline)
	  -fl , --file_loc      Directory in which the data are stored, if set to None, then use the
	                        default data directory ~/.data/geneplexus (default: None)
	  -n , --net_type       Network to use. The choices are: {BioGRID, STRING, IMP} (default: STRING)
	  -f , --features       Types of feature to use. The choices are: {SixSpeciesN2V} (default:
	                        SixSpeciesN2V)
	  -st , --sp_trn        Species of training data The choices are: {Human, Mouse, Fly, Worm,
	                        Zebrafish, Yeast} (default: Human)
	  -sr , --sp_res        Species of results data The choices are: {Human, Mouse, Fly, Worm,
	                        Zebrafish, Yeast}. If more than one species make comma seaprated.
	                        (default: Human)
	  -gt , --gsc_trn       Geneset collection used to generate negatives. The choices are: {GO,
	                        Monarch, Mondo, Combined} (default: Combined)
	  -gr , --gsc_res       Geneset collection used for model similarities. The choices are: {GO,
	                        Monarch, Mondo, Combined}. If more than one gsc can be comma spearated.
	                        (default: Combined)
	  -in , --input_negatives
	                        Input negative gene list (.txt) file. (default: None)
	  -l , --log_level      Logging level. The choices are: {CRITICAL, ERROR, WARNING, INFO, DEBUG}.
	                        Set to CRITICAL for quietest logging. (default: INFO)
	  -ad, --auto_download  When added turns on autodownloader which is off by default. (default:
	                        False)
	  --clear-data          When added will allow user to interactively clear file_loc data and exit.
	                        (default: False)
	  --do_clustering       When added cluster_input() function will be run. (default: False)
	  --skip-mdl-sim        When added make_sim_dfs() will not be run (default: False)
	  --skip-sm-edgelist    When added make_small_edgelist() will not be run (default: False)
	  -cm , --clust_method
	                        Sets the clustering method in cluster_input(). The choices are: {louvain,
	                        domino} (default: louvain)
	  -cmin , --clust_min_size
	                        Sets the minimum size of clusters allowed in cluster_input(). (default:
	                        15)
	  -cw, --clust_weighted
	                        When added will set clust_weight argument to False in cluster_input().
	                        (default: True)
	  -ck , --clust_kwargs
	                        Sets the clustering keyword arguments in cluster_input(). (default:
	                        {'louvain_max_size': 70, 'louvain_max_tries': 3, 'louvain_res': 1,
	                        'louvain_seed': 123, 'domino_res': 1, 'domino_slice_thresh': 0.3,
	                        'domino_n_steps': 20, 'domino_module_threshold': 0.05, 'domino_seed':
	                        123})
	  -lk , --logreg_kwargs
	                        Set the logistic regression keyword arguments in fit(). (default:
	                        {'max_iter': 10000, 'solver': 'lbfgs', 'penalty': 'l2', 'C': 1.0})
	  -s, --scale           When added, will set scale to True in fit(). See docs for more info of
	                        when this is good to do. (default: False)
	  -mnp , --min_num_pos
	                        Minimum umber of genes needed to fit a model in fit(). (default: 15)
	  -mnpcv , --min_num_pos_cv
	                        Minumum number of genes needed to do cross validation in fit(). (default:
	                        15)
	  -nf , --num_folds     Number of folds to do for cross validation in fit(). (default: 3)
	  -nv , --null_val      Value to use when CV can't be done in fit(). (default: None)
	  -rs , --random_state
	                        Random state value to use in fit(). (default: 0)
	  -cv, --cross_validate
	                        When added, will set cross validate to False in fit(). (default: True)
	  -nn , --num_nodes     Number of nodes in make_small_edgelist(). (default: 50)
	  -od , --output_dir    Output directory with respect to the repo root directory used in
	                        save_class(). if set to None, then use the default output directory
	                        ~/.data/geneplexus_outputs/results (default: None)
	  -svt , --save_type    Which file saving method to use in save_class(). The choices are: {all,
	                        results_only} (default: all)
	  -z, --zip-output      When added, zip_ouput is set to True in save_class(). (default: False)
	  -o, --overwrite       When added, overwrite is set to True in save_class(). (default: False)


The output file structure is as follows. This is for `--save_type all`, if `--save_type results_only` is
used then only select files will be saved.

* ``my_result/`` Output directory

   * ``geneplexus.log`` The logger messages.
   * ``top_level_config.json`` Contains configuration infomration for GenePlexus class.
   * ``df_convert_out.tsv`` Table showing conversion of input genes to Entrez IDs for all networks. (see :meth:`geneplexus.GenePlexus.load_genes`)
   * ``Model Directories`` Folders containing information for each of the trained models. `All-Genes` for full input gene list and `Cluster-N` for each cluster if clustering was performed.

      * ``clf.joblib`` Serialized version of the trained model.
      * ``std_scale.joblib`` Serialized version of the standard scaler used (``None`` if `scale=False` in :meth:`geneplexus.GenePlexus.fit`).
      * ``model_level_config.json`` Contains configuration information specific to each model including evaluation metrics and positive, megative and neutral genes, and model weights.
      * ``df_convert_out_for_model.tsv`` Table showing conversion of input genes for each model. (see :meth:`geneplexus.GenePlexus.fit`)
      * ``Result Directories`` Folders containing results for each ``sp_res`` and ``gsc_res`` combination

	     * ``df_probs.tsv`` Top predicted genes related to the input gene list. (see :meth:`geneplexus.GenePlexus.predict`)
	     * ``df_sim.tsv`` Similarity of model trained on user gene list to models trained on known gene sets. (see :meth:`geneplexus.GenePlexus.make_sim_dfs`)
	     * ``df_edge.tsv`` Edgelist (Entrez ID) of subgraph induced by top predicted genes. (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
	     * ``df_edge_sym.tsv`` Edgelist (Symbol) of subgraph induced by top predicted genes. (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
	     * ``isoloated_genes.txt`` List of top predicted genes (Entrez ID) that have no edges in the network. (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
	     * ``isoloated_genes_sym.txt`` List of top predicted genes (Symbol) that have no edges in the network. (see :meth:`geneplexus.GenePlexus.make_small_edgelist`)
