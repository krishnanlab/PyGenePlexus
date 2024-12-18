PyGenePlexus FAQs
=====================

Frequently Asked Questions
--------------------------

- `How are positive and negative genes determined? <#q1>`_
- `How to only use negatives that I provide? <#q2>`_
- `Why is there a minimum number of genes I need? <#q3>`_
- `How can I get results from the first version of the web server? <#q4>`_

How are positive and negative genes determined?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. _q1:

In the supervised machine learning model, any gene from the user-supplied
gene list that is able to be converted to an Entrez ID and is also in the network is
considered part of the positive class.

Genes are assigned to the negative class based on the chosen Geneset Context. The default Geneset
Context is Combined, which used all available geneset collections.

GenePlexus then automatically selects the genes in the negative class by:

#. Considering the total pool of possible negative genes to be any gene that has an annotation to at least one of the terms in the selected geneset collection.
#. Retaining all terms in the selected geneset collection that have between 10 and 200 genes annotated to them.
#. Removing genes that are in the positive class.
#. Performing a hypergeometric test between the genes in the positive class and the lists of genes annotated to every term in the selected geneset collection. If the value of this hypergeometric test is less than 0.05, all genes from the given term are also removed from the pool of possible negative genes.
#. Declaring all the remaining genes in the pool of possible negative genes as the negative class.

How to only use negatives that I provide?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. _q2:

If a user wishes to circumvent GenePlexus from automatically finding negative genes, and only use negatives supplied by the user, this can only be done using the API method. An example is provided below. Note that this is highly not recommended to do.

.. code-block:: python

	import geneplexus
	import numpy as np

	gsc_trn = "Combined"
	gsc_res = "Combined"
	features = "SixSpeciesN2V"
	net_type = "STRING"
	sp_trn = "Human"
	sp_res = "Human"

	input_genes = ['10102', '10131', '10140', '10209', '10247', '10394', '10492', '10600', '10605']
	input_negatives = ['92482', '582', '583', '585', '129880', '55212', '27241', '123016']

	gp = geneplexus.GenePlexus(file_loc = fp_data,
	                           gsc_trn = gsc_trn,
	                           gsc_res = gsc_res,
	                           features = features,
	                           net_type = net_type,
	                           sp_trn = sp_trn,
	                           sp_res = sp_res,
	                           log_level="INFO")

	gp.load_genes(input_genes)
	gp.load_negatives(input_negatives)
	myneg = gp.convert_ids_negatives

	gp.pos_genes_in_net, gp.genes_not_in_net, gp.net_genes = geneplexus._geneplexus._get_genes_in_network(
	    gp.file_loc,
	    gp.sp_trn,
	    gp.net_type,
	    gp.convert_ids,
	)

	gp.negative_genes = np.intersect1d(myneg,gp.net_genes)
	gp.neutral_gene_info = {}

	gp.mdl_weights, gp.probs, gp.avgps = geneplexus._geneplexus._run_sl(
	    gp.file_loc,
	    gp.sp_trn,
	    gp.sp_res,
	    gp.net_type,
	    gp.features,
	    gp.pos_genes_in_net,
	    gp.negative_genes,
	    gp.net_genes,
	    logreg_kwargs=None,
	    min_num_pos_cv=15,
	    num_folds=3,
	    null_val=None,
	    random_state=0,
	    cross_validate=True,
	)

	gp.df_probs = geneplexus._geneplexus._make_prob_df(
	    gp.file_loc,
	    gp.sp_trn,
	    gp.sp_res,
	    gp.net_type,
	    gp.probs,
	    gp.pos_genes_in_net,
	    gp.negative_genes,
	)

	gp.make_sim_dfs()

	gp.make_small_edgelist()

Why is there a minimum number of genes I need?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. _q3:

GenePlexus was evaluated on gene sets with at least 10 positive training genes and at least 10 positive testing genes. GenePlexus did not display a significant variation in performance based on the number of genes when the number of training positive genes was small, therefore this default has been set to five to allow more models to be trained.

How can I get results from the first version of the web server?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. _q4:

A web server has been built on top of this python package. If you were utilizing the original version of the web server which gave you a link to your job results, and you would like to retrieve those results, please use the contact info at the bottom of the `Krishnan Lab <https://www.thekrishnanlab.org>`_ web site to let us know. We will help you retrieve these results. Alternatively, `v.1.0.1 <https://pygeneplexus.readthedocs.io/en/v1.0.1/>`_ of the python package can be used to reproduce those results.
