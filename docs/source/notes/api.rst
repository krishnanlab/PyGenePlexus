PyGenePlexus API
================

Download datasets
-----------------

Download necessary files to directory ``data/`` for all tasks for network
[BioGRID]_, using :term:`Embedding` as features, and the geneset collections
(:term:`GSC`\s) [GO]_ and [DisGeNet]_.

.. code-block:: python

   import geneplexus
   geneplexus.download.download_select_data("data", tasks="All", networks="BioGRID",
                                            features="Embedding", GSCs=["GO", "DisGeNet"])

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

Run the PyGenePlexus pipeline
-----------------------------
