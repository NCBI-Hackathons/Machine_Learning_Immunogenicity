# Machine_Learning_Immunogenicity

This is a repo for the Machine Learning Immunogenicity Team in the August 2016 NCBI Hackathon

This project looks into the application of Machine Learning (ML) techniques in the prediction of Immunogenicity (Categorical; Positive or Negative) based on Peptide and its associated Amino Acid properties. The Immune Epitode Database (IEDB) peptide data has been used as part of this study. Additionally, R package "Peptides" has been used to compute the Amino Acid properties and mashup with Peptide data to enable the use of ML algorithms for Immunogenicity analysis - particularly, the algorithms that are more efficient with numeric and categorical data instead of string sequence.

Tensorflow is an open source software library ML that provides linear regression and classification algorithms (open sourced by Google in Nov 2015) for multi-dimensional arrays (aka “Tensors”). K-fold cross-validation as well as hold-out of test data was used to train and test the generated models.

Initial application of Logistic Regression (LR) and Neural Networks (NN) looks promising with approximately 82% and 90% predictive accuracy respectively. Note: Further cross-validation and rigorous analysis needs to be performed to validate these performance metrics. Various other ML algorithms such as variants of Neural Networks such as Convoluted NN, RESNET, MUST-CNN as well as Random Forest, Bayesian Networks should be considered as part of future work. 

The following are provided:
* R scripts for data wrangling of IEDB data and mashup with Amino Acid properties
* Python Notebook for application of Logistic Regression and Neural Networks using Tensorflow

As part of initial results, the convergence of predictive accuracy for Neural Network is presented below.
![alt tag](https://github.com/NCBI-Hackathons/Machine_Learning_Immunogenicity/blob/master/pics/PredictiveAccuracy_NN_InitialFindings.PNG)
