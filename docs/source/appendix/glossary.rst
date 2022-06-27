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

   Adjacency
      Adjacency matrix representation of a network. Let :math:`G = (V, E, w)`
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
      the graph :math:`G` with some propability :math:`\alpha \in (0, 1)` to
      restart, i.e., teleporting back to the starting node.

      More specifically, let :math:`W_D = A D^{-1}` be the random walk matrix
      (column normalized), where :math:`D` is a diagonal matrix of node
      degrees: :math:`D_{i,i} = \text{deg}(i) = \sum_{j \in V} A_{i,j}`.
      Furthermore, let :math:`y \in \mathbb{R}^{|V|}` be a probability
      distribution of initial "heat" in each node. Then, the one hop random
      walk (or propagation) is :math:`\text{PROP}(G, y) = W_D y`.

      Finally, we can iteratively compute the random walk (or heat)
      distribution :math:`y^{(t+1)}` at :math:`t+1` step as

      .. math::

         y^{(t+1)} = \alpha y^{(0)} + (1 - \alpha) \text{PROP}(G, y^{(t)})

      And the RWR distribution is taken as
      :math:`\hat y = \lim_{t \to \infty} y^{(t)}`

   Influence
      Influence matrix representation of a network. It is the closed form
      solution of the :term:`RWR`. In particular, the influence matrix is

      .. math::

         F = \alpha (I - (1 - \alpha)W_D)^{-1}

      Then, given any initial heat distribution :math:`y^{(0)}`, the solution
      to the RWR is

      .. math::

         \hat y = \lim_{t \to \infty} y^{(t)} = F y^{(0)}

   Embedding
      Network embeddings generated using [node2vec]_.

   GSC
      A gene set collection (GSC) is a set of gene sets, each of which defines
      a set of positive genes for a specific term (i.e., a label). Currently
      supported GSCs are [GO]_ and [DisGeNet]_.
