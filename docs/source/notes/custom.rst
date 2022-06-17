Using custom networks
=====================

Required data files
-------------------

PyGenePlexus provides helper functions for setting up the necessary files in
order to run the pipeline using a custom network. The required files are:

#. ``NodeOrder_{your_net_name}.txt`` Network node ordering
    A text file with a single column containing all genes present in the
    network. The ordering of nodes in this file serves as the index map for the
    network feature data.
#. ``Data_{feature_type}_{your_net_name}.npy`` A numpy array of the chosen network
    representation (rows are genes ordered by NodeOrder file, columns are features).
#. ``GSC_{gsc_type}_{your_net_name}_GoodSets.json`` Filtered GSC for the network
    A subsetted :term:`GSC` where only the genes present in the network are
    considered. After the intersection, any gene set with size larger than
    ``max_size`` or smaller than ``min_size`` is discarded.
#. ``GSC_{gsc_type}_{your_net_name}_universe.txt`` Common genes in GSC and the network
    All the genes that are present both in the network and in the :term:`GSC`.

Set up files using :mod:`geneplexus.custom`
-------------------------------------------------

Suppose you have an :term:`edgelist` file located at ``path/to/net.edg``
that represents your network (named ``your_net_name``) of interest. And
your working data files are located at ``path/to/data/``. You need to first
set up the working files by

.. code-block:: python

   from geneplexus import custom

   custom.edgelist_to_nodeorder("path/to/net.edg", "path/to/data",
                           "your_net_name")

   # Create adjacency feature representation
   custom.edgelist_to_matrix("path/to/net.edg", "path/to/data",
                             "your_net_name", "Adjacency")

   # Create influence feature representation
   custom.edgelist_to_matrix("path/to/net.edg", "path/to/data",
                             "your_net_name", "Influence", beta=0.85)

   # Set up GO GSC with a minimum gene set size of five
   custom.subset_gsc_to_network("path/to/data", "your_net_name", "GO",
                                min_size=5)

.. Note::

   The custom network setup above only needs to be done once.

After the necessary files are generated, one can then run the GenePlexus
pipeline (see :ref:`PyGenePlexus API`) using the custom network

.. code-block:: python

   from geneplexus import GenePlexus

   gp = GenePlexus("path/to/data/", "your_net_name", "Adjacency", "GO")
   ...
