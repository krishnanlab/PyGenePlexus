Determing Postitive and Negative Samples
========================================

PyGenePlexus uses a supervised machine learning model for predicting the
association of all the genes in the network to the user supplied gene set. To
build the classification boundary the model requires both positive and negative
training examples. The positive set of genes is any gene from the user-supplied
gene list that is able to be converted to an Entrez ID and found in the chosen
network. The user can then choose if they want to define genes in the negative
class based on one of two gene set collections, biological processes from the
[GO]_ or diseases from [DisGeNet]_, based on whether the input genes better
represent a cellular process/pathway or a disease. GenePlexus then
automatically selects the genes in the negative class by:

#. Consider the total pool of possible negative genes to be any gene that has
   an annotation to at least one of the terms in the selected gene set
   collection
#. Remove genes that are in the positive class.
#. For every term in a gene set collection, we perform a one-sided Fisher’s
   exact test between the genes in the positive class and the genes annotated
   to the given term. If the p-value of the test is less than 0.05, all genes
   from the given term are also removed from the pool of possible negative
   genes.
#. The remaining genes in the pool of possible negative genes are used in the
   negative class. Note that most genes in the network are not contained in the
   positive class or negative class and are considered as part of the unlabeled
   class.
