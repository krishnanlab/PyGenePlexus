PyGenePlexus with R
===================

PyGenePlexus can be run in R using the ``reticulate`` pacakge.

.. code-block:: r

   # Load the libraries
   library(reticulate)
   # Path to python environment with PyGenePlexus installed
   use_python("PATH/TO/PYTHON")
   geneplexus <- import("geneplexus")
   # Note: Below TRUE needs to be in R style boolean
   gp <- geneplexus$GenePlexus()
   input_genes <- c("ARL6", "BBS1", "BBS10", "BBS12", "BBS2", "BBS4",
	            "BBS5", "BBS7", "BBS9", "CCDC28B", "CEP290", "KIF7",
		    "MKKS", "MKS1", "TRIM32", "TTC8", "WDPCP")
   # Note: Below need to use $ and not . to access pacakge methods
   gp$load_genes(input_genes)
   gp$fit()
   gp$predict()
   df_probs <- gp$model_info[["All-Genes"]]$results[["Human-Combined"]]$df_probs

   write.csv(df_probs,file="probs.csv")