# Machine_Learning_Immunogenicity

This is a repo for the Machine Learning Immunogenicity Team in the August 2016 NCBI Hackathon

This project looks into the application of Machine Learning (ML) techniques in the prediction of Immunogenicity (positive or negative) based on Peptide and its associated Amino Acid properties. The Immune Epitode Database (IEDB) peptide data has been used as part of this study. R package Peptides has been used to compute the amino acid properties.

Tensorflow is an open source software library ML that provides linear regression and classification algorithms (open sourced by Google in Nov 2015) for multi-dimensional arrays (aka “Tensors”). K-fold validation and hold-out test data was used to test the generated models.

Initial application of Logistic Regression (LR) looks promising with approximately 82% accuracy. Various other ML algorithms such as Neural Networks (NN) and its variant Convoluted NN, RESNET, MUST-CNN as well as Random Forest, Bayesian Networks will be evaluated during the course of the effort.
