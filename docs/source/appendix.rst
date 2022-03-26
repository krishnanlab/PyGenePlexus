Appendix
========

Glossary
--------

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

   Adjacency
      Adjacency matrix representation of a network. Let :math:`G = (V, E, s)`
      be a connected weighted undirected graph with node index set :math:`V`,
      edge set :math:`E`, and an edge weight mapping :math:`w`. The adjacency
      matrix is

      .. math::

         A_{i,j} = \begin{cases}
            w(i, j) & \text{if } (i, j) \in E\\
            0 & \text{if } (i, j) \notin E\\
         \end{cases}

   RWR
      Random walk with restart (RWR) is the process of iteratively walking on
      the graph :math:`G` with some propability :math:`\beta \in (0, 1)` to
      restart, i.e., teleporting back to the starting node.

      More specifically, let :math:`P = A D^{-1}` be the random walk matrix
      (column normalized), where :math:`D` is a diagonal matrix of node
      degrees: :math:`D_{i,i} = \text{deg}(i) = \sum_{j \in V} A_{i,j}`.
      Furthermore, let :math:`y \in \mathbb{R}^{|V|}` be a probability
      distribution of initial "heat" in each node. Then, the one hop random
      walk (or propagation) is :math:`\text{PROP}(G, y) = P y`.

      Finally, we can iteratively compute the random walk (or heat)
      distribution :math:`y^{(t+1)}` at :math:`t+1` step as

      .. math::

         y^{(t+1)} = \beta y^{(0)} + (1 - \beta) \text{PROP}(G, y^{(t)})

      And the RWR distribution is taken as
      :math:`\hat y = \lim_{t \to \infty} y^{(t)}`

   Influence
      Influence matrix representation of a network. It is the closed form
      solution of the :term:`RWR`. In particular, the influence matrix is

      .. math::

         F = \beta (I - (1 - \beta)P)^{-1}

      Then, given any initial heat distribution :math:`y^{(0)}`, the solution
      to the RWR is

      .. math::

         \hat y = \lim_{t \to \infty} y^{(t)} = F y^{(0)}

   Embedding
      Node2vec embeddings generated using [node2vec]_.

   GSC
      A gene set collection (GSC) is a set of gene sets, each of which defines
      a set of positive genes for a specific term (i.e., a label). Currently
      supported GSCs are [GO]_ and [DisGeNet]_.

References
----------

.. [GenePlexus]
   * Liu, R., Mancuso, C.A. *et al.* (2020) `Supervised learning is an accurate method for network-based gene classification <https://doi.org/10.1093/bioinformatics/btaa150>`_. *Bioinformatics* 36, 3457–3465
.. [GO]
   * Ashburner, M., Ball, C., Blake, J. *et al.* (2000) `Gene Ontology: tool for the unification of biology <https://doi.org/10.1038/75556>`_. *Nat Genet* 25, 25–29
   * The Gene Ontology Consortium (2021) `The Gene Ontology resource: enriching a GOld mine <https://doi.org/10.1093/nar/gkaa1113>`_. *Nucl. Acids Res.* 49, D325-D334
.. [DisGeNet]
   * Piñero, J. *et al.* (2019) `The DisGeNET knowledge platform for disease genomics: 2019 update <https://doi.org/10.1093/nar/gkz1021>`_. *Nucl. Acids Res.* 48, D845-D855
   * Piñero, J. *et al.* (2017) `DisGeNET: a comprehensive platform integrating information on human disease-associated genes and variants <https://doi.org/10.1093/nar/gkw943>`_. *Nucl. Acids Res.* 45, D833-D839
   *  Piñero, J. *et al.* (2015) `DisGeNET: a discovery platform for the dynamical exploration of human diseases and their genes <https://doi.org/10.1093/database/bav028>`_. *Database*
.. [BioGRID]
   * Oughtred, R. *et al.* (2020) `The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions <https://doi.org/10.1002/pro.3978>`_. *Protein Sci.*
   * Stark, C. *et al.* (2006) `BioGRID: a general repository for interaction datasets <https://doi.org/10.1093/nar/gkj109>`_. *Nucl. Acids Res.* 34, D535–D539
.. [STRING]
   * Szklarczyk, D. *et al.* (2015) `STRING v10: protein–protein interaction networks, integrated over the tree of life <https://doi.org/10.1093/nar/gku1003>`_. *Nucl. Acids Res.* 43, D447–D452
.. [STRING-EXP]
   * Szklarczyk, D. *et al.* (2015) `STRING v10: protein–protein interaction networks, integrated over the tree of life <https://doi.org/10.1093/nar/gku1003>`_. *Nucl. Acids Res.* 43, D447–D452
.. [GIANT-TN]
   * Greene, C., Krishnan, A., Wong, A. *et al.* (2015) `Understanding multicellular function and disease with human tissue-specific networks <https://doi.org/10.1038/ng.3259>`_. *Nat Genet* 47, 569–576
.. [node2vec]
   * Grover, A., Leskovec, J. (2016) `Node2vec: Scalable Feature Learning for Networks <https://doi.org/10.1145/2939672.2939754>`_. *KDD '16* 855–864
