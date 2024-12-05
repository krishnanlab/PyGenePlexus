PyGenePlexus FAQs
=====================

Frequently Asked Questions
--------------------------

**How are positive and negative genes determined?**

In the supervised machine learning model, any gene from the user-supplied
gene list that is able to be converted to an Entrez ID and is also in the network is
considered part of the positive class.

Genes in the negative class based on the chosen Geneset Context. The default Geneset
Context is Combined, which used all avilable geneset collections.

GenePlexus then automatically selects the genes in the negative class by:

#. Considering the total pool of possible negative genes to be any gene that has an annotation to at least one of the terms in the selected geneset collection.
#. Retaining all terms in the selected geneset collection that have between 10 and 200 genes annotated to them.
#. Removing genes that are in the positive class.
#. Performing a hypergeometric test between the genes in the positive class and the lists of genes annotated to every term in the selected geneset collection. If the value of this hypergeometric test is less than 0.05, all genes from the given term are also removed from the pool of possible negative genes.
#. Declaring all the remaining genes in the pool of possible negative genes as the negative class.