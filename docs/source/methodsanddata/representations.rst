Network Representations
=======================

We utilize three distinct representations of molecular networks: the adjacency matrix, an influence matrix, and low-dimensional node embeddings. Let :math:`G = (V,E,W)` denote an undirected molecular network, where :math:`V` is the set of vertices (genes), :math:`E` is the set of edges (associations between genes), and :math:`W` is the set of edge weights (the strengths of the associations). :math:`G` can be represented as a weighted adjacency matrix :math:`A_{i,j}=W_{i,j}`, where :math:`A=R^{VxV}`. :math:`G` can also be represented as an influence matrix, :math:`F=R^{VxV}`, which can capture both local and global structure of the network. :math:`F` was obtained using a random walk with restart transformation kernel, F=a[I-(1-a)W_{D}]^{-1}

where,  is the restart parameter, I is the identity matrix, and WD is the degree weighted adjacency matrix given by WD=AD-1, where DRVV is a diagonal matrix of node degrees. A restart parameter of 0.85 was used for every network.

G can also be transformed into a low-dimensional representation through the process of node embedding. In this study we used the node2vec algorithm, which borrows ideas from the word2vec algorithm from natural language processing. The objective of node2vec is to find a low-dimensional representation of the adjacency matrix, ERVd, where dV. This is done by optimizing the following log-probability objective function: E=uVlogPrNsu|eu

where Nsu is the network neighborhood of node u generated through a sampling strategy S, and euRd is the feature vector of node u. In node2vec, the sampling strategy is based on random walks that are controlled using two parameters p and q, in which a high value of q keeps the walk local (a breadth-first search), and a high value of p encourages outward exploration (a depth-first search). The values of p and q  were both set to 0.1 for every network.


