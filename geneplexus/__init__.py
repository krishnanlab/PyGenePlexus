"""A Python package for GenePlexus.

[GenePlexus]_ is a network based supervised learning method for gene
classification. Given a list of input genes and a selection of background gene
set collection (:term:`GSC`), it trains a logistic regression model using one
of three network derived features (:term:`adjacency`, :term:`influence`, or
:term:`embedding`) and generate genome wide predictions about genes that are
functionally similar to the input gene list.

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
