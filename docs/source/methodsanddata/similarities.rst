Calulating Model Similarities
=============================

A unique feature of PyGenePlexus is proving some interpretation of the machine learning model trained on the user supplied gene set. This is done by comparing the weights of that trained model to the weights from thousands of other models pretrained on known gene sets of biological processes from [GO]_ and diseases from [DisGeNet]_.

Section 6.1: Pre-training models
--------------------------------

The first step in this process is to train models for each known gene set. For each gene set in either the Gene Ontology or DisGenet collection a model is trained for every combination of network (BioGRID, STRING, STRING-EXP, GIANT-TN), feature type (Adjacency, Influence, Embedding) and way of selecting negatives (GO, DisGeNet) and the weights of these trained models are saved.

The next step is building up matrices that will be used for doing background correction of the final similarities presented on the web-server. To accomplish this, we generate a correction matrix, CRNT, where N is number of terms in the gene set collection that the user chose to build the negative class and T is the number of terms in the gene set collection of the target table. For example, if the user chose to build the negative class based on biological processes in GO and the output result table on the web-server is displaying similarity of the user trained model to diseases in DisGeNet, then the rows of C would correspond to biological process terms and the columns would correspond to disease terms. An element in the correction matrix is given by Ci,j=S(wNi,wTj), where w is the vector of weights from a trained model and S is a function that captures similarity between the two weight vectors. In this work, we use the cosine similarly as our similarity metric. A separate correction matrix is generated for all combinations of network, feature type, negative selection method and target table. We note that this requires training over >10,000 machine learning models where a model can have up to 25,689 weights and use thousands of training examples.

Section 6.2: Getting similarities to user trained model
-------------------------------------------------------

After the user submits a job, a custom machine learning model is trained. Once trained, the GenePlexus web-server computes the cosine similarity of the weights from the user model to the weights of each term in the target table, where qR1T is given by qj=S(wU,wTj) and wU are the weights of the user model. This vector q is then appended as the last row of the corresponding correction matrix, CN+1, j=qj. This is done separately for each target table.

Section 6.3: Background correction
----------------------------------

The background correction is done in two parts. First, a z-score is calculated across all scores for the user genes, which is given by,

                                                                     zqj=max0,CN+1, j-N+1N+1,

Where N+1 and N+1 are the mean and standard deviation calculated across the N+1 row of C, respectively. Additionally, a z-score is calculated to correct for any bias in the negative gene selection, which is given by,

                                                                      zTj=max0,CN+1, j-jj,

Where j and j are the mean and standard deviation calculated across the jth column of C, respectively. The final scores presented on the GenePlexus web-server are the l2-norm of the above z-scores given by,

                                                                          zj=zqj2+zTj2.



