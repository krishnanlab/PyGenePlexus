ML Model
========

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
