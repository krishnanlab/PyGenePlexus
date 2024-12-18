"""PyGenePlexus is a Python package for running the [GenePlexus]_ model.

.. currentmodule:: geneplexus.GenePlexus

PyGenePlexus enables researchers to predict genes similar to an uploaded
geneset of interest based on patterns of connectivity in genome-scale
molecular interaction networks, with the ability to translate these
findings across species.

.. figure:: ../figures/mainfigure.png
  :scale: 20 %
  :align: center
  :alt: My Text

  Overview of PyGenePlexus

Given a list of input genes and a geneset collection (:term:`GSC`) to help
select negative examples, the package trains a logistic regression model
using node embeddings as features and generates the following outputs,
either in the same species as the input genes or translated to a model
species.

#. Genome-wide prediction of how **functionally similar** a gene is to the \
input gene list. Evaluation of the model is provided by performing \
k-fold cross validation. The default is 3-fold cross validation when a \
minimum of 15 input genes are supplied. These parameters can be changed when accessing \
the Python class. PyGenePlexus does not enforce a minimum or maximum number \
of genes, and we note evaluations of the model were carried out for gene sets \
ranging between 5 and 500 genes. See :meth:`fit_and_predict`
#. (Optional) Interpretability of the model is provided by comparing the \
model trained on the user gene set to models pretrained on 1000's of known \
gene sets from [GO]_ bioloigcal proceses, [Monarch]_ phenotypes and [Mondo]_ diseases. See \
:meth:`make_sim_dfs`
#. (Optional) Interpretability of the top predicted genes is provided by \
returning their network connectivity. :meth:`make_small_edgelist`

.. note::

    **Links to other GenePlexus products**

    * `Cross Species GenePlexus Paper (GenePlexusZoo) <https://doi.org/10.1371/journal.pcbi.1011773>`_
    * `GenePlexus Web Server <https://www.geneplexus.net>`_
    * `GenePlexus Web Server Paper \
<https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac335/6586869?login=true>`_
    * `GenePlexus Benchmarking Paper \
<https://academic.oup.com/bioinformatics/article/36/11/3457/5780279>`_
    * `GenePlexus Benchmarking Paper Repo \
<https://github.com/krishnanlab/GenePlexus>`_
    * `PyGenePlexus PyPi Package <https://pypi.org/project/geneplexus/>`_
    * `PyGenePlexus GitHub Repo <https://github.com/krishnanlab/PyGenePlexus>`_
    * `PyGenePlexus Paper \
<https://www.biorxiv.org/content/10.1101/2022.07.02.498552v1.abstract>`_
    * `GenePlexus Web Server Repo <https://github.com/krishnanlab/geneplexus-app-v2>`_

Quick start
-----------

PyGenePlexus comes with an easy to use command line interface (CLI) to run the
full GenePlexus pipeline given an input gene list. Go get started, install via
pip and run a quick example as follows.

.. code-block:: bash

    pip install geneplexus
    geneplexus --input_file my_gene_list.txt --output_dir my_result

Note that you need to supply the ``my_gene_list.txt`` file, which is a line
separated gene list text file  (NCBI Entrez IDs, Symbol or Ensembl IDs are
accepted). An example can be found on the
`GitHub page <https://github.com/krishnanlab/PyGenePlexus>`_ under
``example/input_genes.txt``. More info can be found in :ref:`PyGenePlexus CLI`.

.. warning::

    All necessary files for a specific selection of parameters (network,
    feature, species, and gene set collection) will be downloaded automatically and
    saved under ``~/.data/geneplexus``. User can also specify the location of
    data to be saved using the ``--output_dir`` argument. **The example
    provided will download files that occupy ~4GB of space.**


Using the API
^^^^^^^^^^^^^

A quick example of generating predictions using an input gene list. More info
can be found in :ref:`PyGenePlexus API`.

.. code-block::

    from geneplexus import GenePlexus
    input_genes = ["ARL6", "BBS1", "BBS10", "BBS12", "BBS2", "BBS4",
                   "BBS5", "BBS7", "BBS9", "CCDC28B", "CEP290", "KIF7",
                   "MKKS", "MKS1", "TRIM32", "TTC8", "WDPCP"]
    gp = GenePlexus(net_type="STRING", features="SixSpeciesN2V",
                    sp_trn="Human", sp_res="Human",
                    gsc_trn="Combined", gsc_res="Combined",
                    input_genes=input_genes, auto_download=True,
                    log_level="INFO")
    df_probs = gp.fit_and_predict()[1]
    print(df_probs.iloc[:10])

.. note::

   v2 of PyGenePlexus is signifcanlty different than v1 and uses
   a different set of backend data, which only includes human data.
   For information of that version see
   https://pygeneplexus.readthedocs.io/en/v1.0.1/


"""
from ._config import config  # noreorder
from . import download
from . import util
from .geneplexus import GenePlexus


__version__ = "2.0.2"
__all__ = ["download", "GenePlexus", "util", "config"]
