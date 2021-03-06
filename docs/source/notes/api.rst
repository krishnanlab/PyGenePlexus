PyGenePlexus API
================

.. currentmodule:: geneplexus.geneplexus

Download datasets
-----------------

Manual download
^^^^^^^^^^^^^^^

The examples below show downloading the data to ``my_data/`` for 1) all tasks for network
[STRING]_, using :term:`Embedding` as features, and the geneset collections (:term:`GSC`\s) [GO]_ and [DisGeNet]_ and 2) the full data.

.. warning::

   **PROCEED WITH CAUTION**
   The first example below (STRING network using Embedding features with GO and
   DisGeNet GSCs) will occupy **~300MB** of space. The second example (full
   download) will occupy **~32GB** of space.

.. code-block:: python

   >>> from geneplexus.download import download_select_data
   >>> download_select_data("my_data", tasks="All", networks="STRING",
   ...                      features="Embedding", gscs=["GO", "DisGeNet"])
   >>> download_select_data("my_data")  # alternatively, download all data at once

See :meth:`geneplexus.download.download_select_data` for more information

**Data options:**

======== =======================================================
Networks [BioGRID]_, [STRING]_, [STRING-EXP]_, [GIANT-TN]_
Features :term:`Adjacency`, :term:`Influence`, :term:`Embedding`
GSCs     [GO]_, [DisGeNet]_
======== =======================================================

.. note::

   The :term:`Influence` and :term:`Adjacency` data representations take the longest time to
   download, from **~10 minutes** up to **an hour** dependeing on the download
   speed. The :term:`Embedding` data representation takes the least amount of time to download
   (within **a minute**).

Auto download
^^^^^^^^^^^^^

Optionally, set the ``auto_download`` key word argument to ``True`` to automatically
download necessary data at initialization of the :class:`GenePlexus` object.

.. code-block:: python

   from geneplexus import GenePlexus
   gp = GenePlexus(net_type="STRING", features="Embedding", gsc="GO", auto_download=True)

.. note::

   The default data location is ``~/.data/geneplexus/``. You can change this by
   setting the ``file_loc`` argument of :class:`GenePlexus`.

Run the PyGenePlexus pipeline
-----------------------------

First, specify the input genes (can have mixed gene ID types, i.e. have any combination of Entrez
IDs, Gene Symbols, or Ensembl IDs).

.. code-block:: python

   input_genes = ["6457", "7037", "3134", "TTC8"," BBS5", "BBS12", ...]

Alternatively, read the gene list from file

.. code-block:: python

   import geneplexus
   input_genes = geneplexus.util.read_gene_list("my_gene_list.txt")

Next, run the pipline using the :class:`GenePlexus` object.

.. code-block:: python

   gp = geneplexus.GenePlexus(net_type="STRING", features="Embedding", gsc="GO")

   # Load input genes and set up positives/negatives for training
   gp.load_genes(input_genes)

   # Train logistic regression model and get genome-wide gene predictions
   mdl_weights, df_probs, avgps = gp.fit_and_predict()

   # Optionally, compute model similarity to models pretrained on GO and DisGeNet gene sets
   df_sim_GO, df_sim_Dis, weights_GO, weights_Dis = gp.make_sim_dfs()

   # Optionally, extract the subgraph induced by the top (50 by default) predicted genes
   df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = gp.make_small_edgelist()

