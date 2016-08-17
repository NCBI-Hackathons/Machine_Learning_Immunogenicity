# Machine_Learning_Immunogenicity

This is a repo for the Machine Learning Immunogenicity Team in the August 2016 NCBI Hackathon

This project looks into the application of Machine Learning (ML) techniques in the prediction of Immunogenicity (Categorical; Positive or Negative) based on Peptide and its associated Amino Acid properties. The Immune Epitode Database (IEDB) peptide data has been used as part of this study. Additionally, R package "Peptides" has been used to compute the Amino Acid properties and mashup with Peptide data to enable the use of ML algorithms for Immunogenicity analysis - particularly, the algorithms that are more efficient with numeric and categorical data instead of string sequence.

Tensorflow is an open source software library ML that provides linear regression and classification algorithms (open sourced by Google in Nov 2015) for multi-dimensional arrays (aka “Tensors”). K-fold cross-validation as well as hold-out of test data was used to train and test the generated models.

Initial application of Logistic Regression (LR) and Neural Networks (NN) looks promising with approximately 82% and 90% predictive accuracy respectively. Various other ML algorithms such as variants of Neural Networks such as Convoluted NN, RESNET, MUST-CNN as well as Random Forest, Bayesian Networks should be considered as part of future work. 

The following are provided:
* R scripts for data wrangling of IEDB data and mashup with Amino Acid properties
* Python Notebook for application of Logistic Regression and Neural Networks using Tensorflow
* Python script for computing binding affinities using several published approaches based on several datasets

## Run predictors on published data

This script allows to run published epitope binding predictors on immunogenicity data from IEDB T-cell and MHC assays (http://www.iedb.org/), as well as data from IMMA2. This script can be used to test whether binding predictions also predict immugenicity (Note: they don't). We use the data interface implemented in the pepdata package from the Hammer Lab (https://github.com/hammerlab/pepdata) and the implementations of various predictors in the Fred 2 framework for computational immunogenomics by Schubert et al. (https://github.com/FRED-2/Fred2 and http://bioinformatics.oxfordjournals.org/content/32/13/2044).

### Usage

    usage: run-data-on-predictor.py [-h] --predictor PREDICTOR --dataset DATASET
                                    [-n N] [--allele ALLELE]
    
    Call epitope predictors on data.
    
    optional arguments:
      -h, --help            show this help message and exit
      -n N                  Number of rows to take from dataset
      --allele ALLELE       Allelle
    
    required arguments:
      --predictor PREDICTOR
                            Epitope predictors [see all with --predictor=list]
      --dataset DATASET     Immunogenic dataset [see all with --dataset=list]


## List available predictors and datasets

    run-data-on-predictor.py  --predictor list --dataset list

Set one of the predictors with --predictor:
['smmpmbec', 'syfpeithi', 'netctlpan', 'smm', 'tepitopepan', 'netmhcii', 'arb', 'pickpocket', 'epidemix', 'unitope', 'netmhciipan', 'comblibsidney', 'netmhcpan', 'calisimm', 'hammer', 'svmhc', 'bimas']

Details from https://bioinformatics.oxfordjournals.org/content/suppl/2016/02/26/btw113.DC1/S1.pdf

| SYFPEITHI | T-cell epitope | (Rammensee, et al., 1999)|
| BIMAS | MHC-I binding | (Parker, et al., 1994)|
| SVMHC | MHC-I binding | (Dönnes and Elofsson, 2002)|
| ARB | MHC-I binding | (Bui, et al., 2005)|
| SMM | MHC-I binding | (Peters and Sette, 2005)|
| SMMPMBEC | MHC-I binding | (Kim, et al., 2009)|
| Epidemix | MHC-I binding | (Feldhahn, et al., 2009)|
| Comblib | MHC-I binding | (Sidney, et al., 2008)|
| PickPocket* | MHC-I binding | (Zhang, et al., 2009)|
| NetMHC* | MHC-I binding | (Lundegaard, et al., 2008)|
| NetMHCpan* | MHC-I binding | (Hoof, et al., 2009)|
| HAMMER | MHC-II binding | (Sturniolo, et al., 1999)|
| TEPITOPEpan | MHC-II binding | (Zhang, et al., 2012)|
| NetMHCII* | MHC-II binding | (Nielsen, et al., 2007)|
| NetMHCIIpan* | MHC-II binding | (Karosiene, et al., 2013)|
| UniTope | T-cell epitope | (Toussaint, et al., 2011)|
| NetCTLpan* | T-cell epitope | (Stranzl, et al., 2010)|

available datasets: iedb.tcell, ideb.mhc, imma2

## Example output


