"""A Python package for GenePlexus.

[GenePlexus]_ is a network based supervised learning method for gene
classification. Given a list of input genes and a selection of background gene
set collection (:term:`GSC`), it trains a logistic regression model using one
of three network derived features (:term:`adjacency`, :term:`influence`, or
:term:`embedding`) and generate the followings

#. Genomewide prediction about genes that are **functionally similar** to the \
input gene list. See :meth:`geneplexus.GenePlexus.fit_and_predict`
#. (Optional) Input gene list simiarltiy with [GO]_ or [DisGeNet]_ based on \
the model coefficients. See :meth:`geneplexus.GenePlexus.make_sim_dfs`
#. (Optional) A subnetwork induced by the top predicted genes. See \
:meth:`make_small_edgelist`

.. code-block::

    >>> from geneplexus import GenePlexus, util
    >>> gp = GenePlexus(net_type="BioGRID", features="Embedding", gsc="GO", \
auto_download=True)
    >>> input_genes = util.read_gene_list("example/input_genes.txt")
    >>> input_genes[:7]
    ["6457", "7037", "57403", "3134", "50807", "93343", "11311"]
    >>> gp.load_genes(input_genes)
    >>> _, df_probs, _ gp.fit_and_predict()
    >>> df_probs.iloc[:4]
       Entrez  Symbol                                              Name\
  Probability Known/Novel Class-Label  Rank
    0   57154  SMURF1       SMAD specific E3 ubiquitin protein ligase 1\
         1.0       Known           P     1
    1    1213    CLTC                              clathrin heavy chain\
         1.0       Known           P     2
    2    4734   NEDD4                 NEDD4 E3 ubiquitin protein ligase\
         1.0       Known           P     3
    3    6714     SRC  SRC proto-oncogene, non-receptor tyrosine kinase\
         1.0       Known           P     4


Currently, GenePlexus come with four networks, including [BioGRID]_ (default),
[STRING]_, [STRING-EXP]_, and [GIANT-TN]_. Prediction using custom network can
also be done, see :ref:`Using custom networks`. However, when using custom
network, the model similarity (against GO and DisGeNet) analysis *cannot* be
done due to the lack to pretrained models.

.. note::

    A webserver for GenePlexus is also available at
    `<https://www.geneplexus.net>`_

"""
from ._config import config  # noreorder
from . import download
from . import util
from . import custom
from .geneplexus import GenePlexus


__all__ = ["download", "GenePlexus", "util", "config", "custom"]
