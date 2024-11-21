PyGenePlexus Data
=====================

Preprocessed Data
-------------------

PyGenePlexus comes with pre-processed data that can be downloaded
using :meth:`geneplexus.download.download_select_data` or directly
from `Zenodo <https://zenodo.org/records/14149956>`_.

**Data options:**

======== =======================================================
Networks [BioGRID]_, [STRING]_, [IMP]_
Species  Human, Mouse, Zebrafish, Worm, Fly, Yeast
Features :term:`SixSpeciesN2V`
GSCs     [GO]_, [Monarch]_, [Mondo]_
======== =======================================================

Custom Data
------------

PyGenePlexus uses a strict naming convention rto read in the files. Users
can supply there own data by generating files in the following formats:

#. ``NodeOrder__{Species}__{Network Name}.txt`` Network node ordering
    A text file with a single column containing all genes present in the
    network for a given species. The ordering of nodes in this file
    serves as the index map for thenetwork feature data.
#. ``Data__{Species}__{Feature Name}__{Network Name}.npy`` Data array
    A numpy array of the chosen network representation (rows are genes
    ordered by NodeOrder file, columns are features).
#. ``GSC__{Species}__{GSC Name}__{Network Name}.json`` Filtered GSC for the network
    A subsetted :term:`GSC` where only the genes present in the network are
    considered.
    ::

       {
         "{Term ID}" # ID of the term : {
            "Name"  : # returns string of term name
            "Genes" : # returns list of genes annotated to term
            "Task"  : # returns ype of GSC the term is from
            }
         "Universe" : # returns list of all genes in GSC
         "Term Order" : # returns list of all term IDs in GSC
       }

#. ``PreTrainedWeights__{Species}__{GSC Name}__{Network Name}__{Feature Name}.json`` Pretrained weights
    Model weights for model trained on terms is selected :term:`GSC`
    ::

       {
         "{Term ID}" # ID of the term : {
            "Name"  : # returns string of term name
            "Weights" : # returns list of model weights
            "PosGenes" : # returns list genes used as positives in the model
            "Task"  : # returns type of GSC the term is from
            }
       }

#. ``Edgelist__{Species}__{Network Name}.edg`` Network edgelist
    :term:`Edgelist` for the given network.
#. ``IDconversion__{Species}__{ID Type}-to-{ID Type}.json`` Gene ID conversions
    These files convert from one gene ID convention to another
    ::

       {
         "{Gene ID}" : # returns list of how Gene ID is converted to other gene ID type
       }

#. ``BioMart__{Species}__{Species}.json`` One-to-one ortholog conversions
    These files convert genes from one species into the corresponding one-to-one ortholog.
    ::

       {
         "{Gene ID}" : # returns string of how gene is converted to its one-to-one ortholog
       }
