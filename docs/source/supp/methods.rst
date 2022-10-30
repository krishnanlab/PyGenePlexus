Methods
=======

ML Model
--------

In GenePlexus, the supervised machine learning model uses the connections of a
user chosen genome-scale molecular network as feature vectors. As described
above, these feature vectors can be one of three representations; Adjacency,
Influence and Embedding. PyGenePlexus uses logistic regression with
l2-regularization as the supervised learning algorithm and is implemented using
the python package
`scikit-learn <https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html>`_.
After training a model using the labeled genes, the trained model is used to
classify all the genes in the chosen network, returning a prediction
probability for these genes that is bounded between 0 and 1. The regularization
parameter is set to 1.0 ny default but can be changed by the user.

Calculating Model Similarities
------------------------------

A unique feature of PyGenePlexus is proving some interpretation of the machine
learning model trained on the user supplied gene set. This is done by comparing
the weights of that trained model to the weights from thousands of other models
pretrained on known gene sets of biological processes from [GO]_ and diseases
from [DisGeNet]_.

Pre-training models
+++++++++++++++++++
The first step in this process is to train models for each known gene set. For
each gene set in either GO or DisGenet collection a model is trained for every
combination of network (BioGRID, STRING, STRING-EXP, GIANT-TN), feature type
(:term:`Adjacency`, :term:`Influence`, :term:`Embedding`) and way of selecting
negatives (GO, DisGeNet) and the weights of these trained models are saved.

The next step is building up matrices that will be used for doing background
correction of the final similarities presented on the web-server. To accomplish
this, we generate a correction matrix, :math:`C{\in}R^{N \times T}`, where
:math:`N` is number of terms in the gene set collection that the user chose to
build the negative class and :math:`T` is the number of terms in the gene set
collection of the target table. For example, if the user chose to build the
negative class based on biological processes in GO and the output result table
on the web-server is displaying similarity of the user trained model to
diseases in DisGeNet, then the rows of :math:`C` would correspond to biological
process terms and the columns would correspond to disease terms. An element in
the correction matrix is given by :math:`C_{i,j}=S((w_N)_i, (w_T)_j)`, where
:math:`w` is the vector of weights from a trained model and :math:`S` is a
function that captures similarity between the two weight vectors. In this work,
we use the cosine similarly as our similarity metric. A separate correction
matrix is generated for all combinations of network, feature type, negative
selection method and target table. We note that this requires training over
>10,000 machine learning models where a model can have up to 25,689 weights and
use thousands of training examples.

Getting similarities to user trained model
++++++++++++++++++++++++++++++++++++++++++
After the user submits a job, a custom machine learning model is trained. Once
trained, the GenePlexus web-server computes the cosine similarity of the
weights from the user model to the weights of each term in the target table,
where :math:`q{\in}R^{1 \times T}` is given by :math:`q_{j}=S(w_U, (w_T)_j)`
and :math:`w_{U}` are the weights of the user model. This vector :math:`q` is
then appended as the last row of the corresponding correction matrix,
:math:`C_{N+1,j}=q_{j}`. This is done separately for each target table.

Background correction
+++++++++++++++++++++
The background correction is done in two parts. First, a z-score is calculated
across all scores for the user genes, which is given by,

.. math::
   z_{q_{j}}=max(0,\frac{C_{N+1,j}-{\mu}(N+1)}{{\sigma}(N+1)}),

Where :math:`{\mu}(N+1)` and :math:`{\sigma}(N+1)` are the mean and standard
deviation calculated across the :math:`N+1` row of :math:`C`, respectively.
Additionally, a z-score is calculated to correct for any bias in the negative
gene selection, which is given by,

.. math::
   z_{T_{j}}=max(0,\frac{C_{N+1,j}-{\mu}(j)}{{\sigma}(j)}),

Where :math:`{\mu}(j)` and :math:`{\sigma}(j)` are the mean and standard
deviation calculated across the :math:`j^{th}` column of :math:`C`,
respectively. The final scores presented on the GenePlexus web-server are the
l2-norm of the above z-scores given by,

.. math::
   z_{j}=\sqrt{z_{q_{j}}^{2}+z_{T_{j}}^{2}}.
