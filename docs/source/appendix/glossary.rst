Glossary
========

.. glossary::
   Edgelist
      An edge-list file (``.edg``) is a representation of the network. It
      contains either two columns (unweighted) or three columns (weighted).
      The first two columns are gene pairs representing edges in the network.
      Since the network is assumed to be undirected, reversed direction is
      unnecessary. When the network is weigthed, the edge weights should be
      encoded in the third column.

      .. code-block::

         gene1  gene2   <weight>

   Embedding
      Network embeddings generated using [node2vec]_.

   GSC
      A gene set collection (GSC) is a set of gene sets, each of which defines
      a set of positive genes for a specific term (i.e., a label). Currently
      supported GSCs are [GO]_, [Monarch]_ and [Mondo]_.

   SixSpeciesN2V
      Features where made by first constucting a six species network, connecting
      genes from the single species networks together using the [EggNOG]_ orthology
      database and then embedding this network using [node2vec]_.
