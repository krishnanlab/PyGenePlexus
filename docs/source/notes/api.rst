PyGenePlexus API
================

.. currentmodule:: geneplexus.geneplexus

Download datasets
-----------------

Download necessary files to directory ``data/`` for all tasks for network
[BioGRID]_, using :term:`Embedding` as features, and the geneset collections
(:term:`GSC`\s) [GO]_ and [DisGeNet]_.

.. code-block:: python

   import geneplexus
   geneplexus.download.download_select_data("data", tasks="All", networks="BioGRID",
                                            features="Embedding", gscs=["GO", "DisGeNet"])

   # Alternatively, to download all data at once
   geneplexus.download.download_select_data("data")

See :meth:`geneplexus.download.download_select_data` for more information

List of data options
    * Networks
        * [BioGRID]_
        * [STRING]_
        * [STRING-EXP]_
        * [GIANT-TN]_
    * Features
        * :term:`Adjacency`
        * :term:`Influence`
        * :term:`Embedding`
    * GSC
        * [GO]_
        * [DisGeNet]_

Auto download
^^^^^^^^^^^^^

Optionally, set to ``auto_download`` key word argument to ``True`` to automatically
download necessary data at initialization of the :class:`GenePlexus` object.

.. code-block:: python

   from geneplexus import GenePlexus
   gp = GenePlexus(net_type="BioGRID", features="Embedding", gsc="GO", auto_download=True)

.. warning::

   :term:`Influence` (followed by :term:`Adjacency`) data takes a long time to
   download, from **~10 minutes** up to **an hour** dependeing on the download
   speed. :term:`Embedding` feature data takes least amount of time to download
   (within **a minute**).

Run the PyGenePlexus pipeline
-----------------------------

First, specify the input genes (can have mixed gene ID types, e.g., Entrez
gene ID, gene Symbol)

.. code-block:: python

   input_genes = ["6457", "7037", "3134", "TTC8"," BBS5", "BBS12", ...]

Alternatively, read the gene list from file

.. code-block:: python

   import geneplexus
   input_genes = geneplexus.util.read_gene_list("my_gene_list.txt")

Next, run the pipline via the :class:`geneplexus.GenePlexus` object. The data
files are stored under the ``~/.data/geneplexus`` directory by default. The
data file location can be changed by setting the ``file_loc`` argument (see
:meth:`geneplexus.GenePlexus.__init__`).

.. code-block:: python

   gp = geneplexus.GenePlexus(net_type="BioGRID", features="Embedding", gsc="GO")

   # Load input genes and set up positives/negatives for training
   gp.load_genes(input_genes)

   # Train logistic regression model and get genomewide gene predictions
   mdl_weights, df_probs, avgps = gp.fit_and_predict()

   # Optionally, compute modle similarity against pretrained models for GO and DisGeNet
   df_sim_GO, df_sim_Dis, weights_GO, weights_Dis = gp.make_sim_dfs()

