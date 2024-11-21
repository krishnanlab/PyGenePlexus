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
   download_select_data("my_data", species = ["Human", "Mouse"]) # download just Human nd Mouse data
   download_select_data("my_data")  # download all data at once

See :meth:`geneplexus.download.download_select_data` for more information

Auto download
^^^^^^^^^^^^^

Optionally, set the ``auto_download`` key word argument to ``True`` to automatically
download necessary data at initialization of the :class:`GenePlexus` object.

.. code-block:: python

   from geneplexus import GenePlexus
   gp = GenePlexus(net_type="STRING", features="SixSpeciesN2V",
                   sp_trn = "Human", sp_res = "Human",
                   auto_download=True)

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

   # Instantiate GenePlexus class with default parameters
   gp = geneplexus.GenePlexus()

   # Load input genes and set up positives/negatives for training
   gp.load_genes(input_genes)

   # Train logistic regression model and get genome-wide gene predictions
   mdl_weights, df_probs, avgps = gp.fit_and_predict()

   # Optionally, compute model similarity to models pretrained on GO and DisGeNet gene sets
   df_sim, weights_dict = gp.make_sim_dfs()

   # Optionally, extract the subgraph induced by the top (50 by default) predicted genes
   df_edge, isolated_genes, df_edge_sym, isolated_genes_sym = gp.make_small_edgelist()

