Construction of the networks
============================

The pre-processed networks are [BioGRID]_, [STRING]_ and [GIANT-TN]_. Detailed information about the network properties and sources can be seen in the table below, with the network construction method and interaction type information coming from [Huang]_ et. al.. BioGRID (version 4.2.191) is a low-throughput network that includes both genetic interactions, as well as physical protein-protein interactions. STRING (version 11.0) is a high-throughput, scored network that aggregates information from many data sources. We used two different STRING networks. First, we used the “combined” network that directly includes database annotations, text-mining, ortholog information, co-expression, and experimental determined interactions (referred to as “STRING”). We also used a subset of edges in STRING that had just the “experiments” data, thus restricting the network to one constructed just from experimental determined interactions in humans (referred to as “STRING-EXP”). For both networks, we used the corresponding relationship scores as edge weights, after scaling them to lie between 0 and 1. The GIANT-TN (version 1.0) network is the tissue-naïve network from GIANT, referred to as the “Global” network on the HumanBase website, and is constructed from both low- and high-throughput data, and includes information from co-expression, non-protein sources, regulatory data, and physical protein-protein interactions. The GIANT-TN network is a fully connected, scored network. To add sparsity to the GIANT-TN network, we removed all edges with scores below 0.01 (equal to the prior in the Bayesian model used to construct the network). We used all edge scores (weights) unless otherwise noted, and the nodes in all networks were mapped into Entrez genes using the [MyGeneInfo]_ database. If the original node ID mapped to multiple Entrez IDs, we added edges between all possible mappings.

.. figure:: ../../figures/networkdata.png
  :scale: 50 %
  :align: center
  :alt: My Text

  Table S1. Information on the molecular networks. LT : low-throughput, HT : high-throughput, G : genetic, P : physical, E : Experimentally determined, DA : database annotations, CE : co-expression, NP : non-protein, R : regulation, TM : text-mining, O : orthologous.


