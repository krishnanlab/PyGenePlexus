Network Representations
=======================

We utilize three distinct representations of molecular networks: the adjacency matrix, an influence matrix, and low-dimensional node embeddings. Let :math:`G = (V,E,W)` denote an undirected molecular network, where :math:`V` is the set of vertices (genes), :math:`E` is the set of edges (associations between genes), and :math:`W` is the set of edge weights (the strengths of the associations). :math:`G` can be represented as a weighted adjacency matrix :math:`A_{i,j}=W_{i,j}`, where :math:`A{\in}R^{V{\times}V}`. :math:`G` can also be represented as an influence matrix, :math:`F{\in}R^{V{\times}V}`, which can capture both local and global structure of the network. :math:`F` was obtained using a random walk with restart transformation kernel,

.. math::
   F=\alpha[I-(1-\alpha)W_{D}]^{-1}

where, :math:`\alpha` is the restart parameter, :math:`I` is the identity matrix, and :math:`W_{D}` is the degree weighted adjacency matrix given by :math:`W_{D}=AD^{-1}`, where :math:`D{\in}R^{V{\times}V}` is a diagonal matrix of node degrees. A restart parameter of 0.85 was used for every network.

:math:`G` can also be transformed into a low-dimensional representation through the process of node embedding. In this study we used the [node2vec]_ algorithm, which borrows ideas from the word2vec algorithm from natural language processing. The objective of *node2vec* is to find a low-dimensional representation of the adjacency matrix, :math:`E{\in}R^{V{\times}d}`, where :math:`d{\ll}V`. This is done by optimizing the following log-probability objective function:

.. math::
   E=\sum_{u{\in}V}{log(Pr(N_{s}(u)|e(u)))}

where :math:`N_{s}(u)` is the network neighborhood of node :math:`u` generated through a sampling strategy :math:`S`, and :math:`e(u){\in}R^{d}` is the feature vector of node :math:`u`. In *node2vec*, the sampling strategy is based on random walks that are controlled using two parameters :math:`p` and :math:`q`, in which a high value of :math:`q` keeps the walk local (a breadth-first search), and a high value of :math:`p` encourages outward exploration (a depth-first search). The values of :math:`p` and :math:`q`  were both set to 0.1 for every network.


