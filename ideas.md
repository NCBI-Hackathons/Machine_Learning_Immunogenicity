* Supervised learning on peptide sequences labeled for immunogenicity
    * Logistic regression - 82%
    * Find another test data set to benchmark against other predictors
    * Submit to IEDB benchmark: http://tools.iedb.org/auto_bench/mhci/weekly/
    * Add features based on peptide physiochemical properties
* Viral sequence data set from Florian's collaborator
* Local benchmarking tool to run our algorithm against others in FRED2
* R/python package for working with IEDB data
* Add mhcflurry to FRED2
* Cancer-associated peptide databases
    * Look for association between peptides and cancer types
    * Add cancer peptides to bacteria/virus peptides and perform unsupervised learning
        * Since cancer peptides are just mutated human peptides, won't the clustering just separate based on species
* Random peptide sequence generator (using RNN?)
