"""A Python package for GenePlexus.

.. currentmodule:: geneplexus.GenePlexus

[GenePlexus]_ is a network based supervised learning method for gene
classification. Given a list of input genes and a selection of background gene
set collection (:term:`GSC`), it trains a logistic regression model using one
of three network derived features (:term:`adjacency`, :term:`influence`, or
:term:`embedding`) and generate the followings

#. Genome-wide prediction about genes that are **functionally similar** to the \
input gene list. See :meth:`fit_and_predict`
#. (Optional) Input gene list simiarltiy with :term:`GSC` based on \
model coefficients. See :meth:`make_sim_dfs`
#. (Optional) A subnetwork induced by the top predicted genes. See \
:meth:`make_small_edgelist`

.. note::

    A webserver for GenePlexus is also available at
    `<https://www.geneplexus.net>`_

Quick start
-----------

PyGenePlexus come with an easy to use command line interface (CLI) to run the
full GenePlexus pipeline given an input gene list. Go get started, install via
pip and run a quick example as follows.

.. code-block:: bash

    pip install geneplexus
    geneplexus -i my_gene_list.txt --output_dir my_result

Note that you need to supply the ``my_gene_list.txt`` file, which is a line
separated gene list (can be Entrez or Symbol) text file. An example can be
found on the `GitHub page <https://github.com/krishnanlab/PyGenePlexus>`_ under
``example/input_genes.txt``. More info can be found in :ref:`PyGenePlexus CLI`

.. warning::

    All necessary files for a specific selection of parameters (network,
    feature, and gene set collection) will be downloaded automatically and
    saved under ``~/.data/geneplexus``. User can also specify the location of
    data to be saved using the ``--output_dir`` argument. **The example
    provided will download files that occupy ~230MB of space.**


Using the API
^^^^^^^^^^^^^

A quick example of generating predictions using an input gene list. More info
can be found in :ref:`PyGenePlexus API`

.. code-block::

    >>> from geneplexus import GenePlexus, util
    >>> input_genes = util.read_gene_list("example/input_genes.txt")
    >>> input_genes[:7]
    ["6457", "7037", "57403", "3134", "50807", "93343", "11311"]
    >>> gp = GenePlexus(net_type="BioGRID", features="Embedding", gsc="GO",
    ...                 input_genes=input_genes, auto_download=True)
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


Supported networks
^^^^^^^^^^^^^^^^^^

Currently, GenePlexus come with four networks, including [BioGRID]_ (default),
[STRING]_, [STRING-EXP]_, and [GIANT-TN]_. Prediction using custom network can
also be done, see :ref:`Using custom networks`. However, when using custom
network, the model similarity (against GO and DisGeNet) analysis *cannot* be
done due to the lack to pretrained models.

"""
from ._config import config  # noreorder
from . import download
from . import util
from . import custom
from .geneplexus import GenePlexus


__all__ = ["download", "GenePlexus", "util", "config", "custom"]
