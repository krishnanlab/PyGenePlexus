PyGenePlexus API
================

.. currentmodule:: geneplexus.geneplexus

Download datasets
-----------------

Manual download
^^^^^^^^^^^^^^^

.. warning::

   **PROCEED WITH CAUTION**
   The first example download below will occupy **~4GB** of space. The second example (full
   download) will occupy **~6.5GB** of space.

.. code-block:: python

   from geneplexus.download import download_select_data
   download_select_data(file_loc=None, species=["Human", "Mouse"]) # download just Human nd Mouse data

See :meth:`geneplexus.download.download_select_data` for more information

Auto download
^^^^^^^^^^^^^

Optionally, set the ``auto_download`` key word argument to ``True`` to automatically
download necessary data at initialization of the :class:`GenePlexus` object.

.. code-block:: python

   from geneplexus import GenePlexus
   gp = GenePlexus(net_type="STRING", features="SixSpeciesN2V",
                   sp_trn="Human", sp_res="Human",
                   file_loc=None, auto_download=True)

.. note::

   The default data location is ``~/.data/geneplexus/``. You can change this by
   setting the ``file_loc`` argument of :class:`GenePlexus`.

Loading an input gene set
-------------------------

First, specify the input genes (can have mixed gene ID types, i.e. have any combination of Entrez
IDs, Gene Symbols, or Ensembl IDs).

.. code-block:: python

   input_genes = ["6457", "7037", "3134", "TTC8"," BBS5", "BBS12", ...]

Alternatively, read the gene list from file

.. code-block:: python

   import geneplexus
   input_genes = geneplexus.util.read_gene_list("my_gene_list.txt")

Example running PyGenePlexus pipeline
-------------------------------------

Next, run the pipline using the :class:`GenePlexus` object.

.. code-block:: python

	import geneplexus
	import json
	import os.path as osp

	# if you downloaded back end data, we can get a set from there to use as input
	# find a set that is large enough to cluster
	# replace path if didn't use deafualt file_loc
	fp_base = osp.expanduser("~")
	fp_full = osp.join(fp_base,
	                   ".data",
	                   "geneplexus",
	                   "PreTrainedWeights__Human__Mondo__STRING__SixSpeciesN2V.json",
	)
	with open (fp_full, "r") as f:
	    disease_gene_sets = json.load(f)
	for aterm in disease_gene_sets:
	    num_genes = len(disease_gene_sets[aterm]["PosGenes"])
	    if num_genes > 100:
	        input_genes = disease_gene_sets[aterm]["PosGenes"]
	        print(f"Disease chosen is {disease_gene_sets[aterm]['Name']}")
	        break

	# initialize GenePlexus
	gp = geneplexus.GenePlexus(file_loc = None,
	                           net_type = "STRING",
	                           sp_trn = "Human",
	                           gsc_trn = "Combined",
	                           sp_res = ["Human", "Mouse"],
	                           gsc_res = ["Combined", "Combined"],
	                           input_genes = input_genes,
	)

	# do clustering and generate all results
	gp.cluster_input()
	gp.fit()
	gp.predict()
	gp.make_sim_dfs()
	gp.make_small_edgelist()

	# get human gene prediction results for full input gene set model
	print(gp.model_info["All-Genes"].results["Human-Combined"].df_probs)
	# get similarties of trainied model to other models trained with human annotations
	print(gp.model_info["All-Genes"].results["Human-Combined"].df_sim)
	# get network connections for the top 50 human genes predcited using full input gene set model
	print(gp.model_info["All-Genes"].results["Human-Combined"].df_edge_sym)

	# get mouse gene prediction results for cluster 1 gene set model
	print(gp.model_info["Cluster-01"].results["Mouse-Combined"].df_probs)
	# get similarties of trainied model to other models trained with mouse annotations
	print(gp.model_info["Cluster-01"].results["Mouse-Combined"].df_sim)
	# get network connections for the top 50 mouse genes predcited using cluster 1 gene set model
	print(gp.model_info["Cluster-01"].results["Mouse-Combined"].df_edge_sym)

	# get log2(auPRC/prior) metric for the full input gene set model
	print(gp.model_info["All-Genes"].avgps)

	# save the class. If output_dir=None will try to save to ~/.data/geneplexus_outputs/results
	gp.save_class(output_dir = None)

- For all items saved in GenePlexus class see :class:`GenePlexus`.
- For structure of save_class output see :ref:`PyGenePlexus CLI <cli>`.