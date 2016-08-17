* Supervised learning on peptide sequences labeled for immunogenicity
    * Logistic regression - 82%
    * Find another test data set to benchmark against other predictors
    * Add features based on peptide physiochemical properties
    * For epitope prediction
        * Add MHC allele to model
        * Submit to IEDB benchmark: http://tools.iedb.org/auto_bench/mhci/join
    * For a given condition (cancer, virus), predict whether an epitope will elicit an immunogenic response
    * How many validated cancer-associated epitopes are predicted from database of known human neo-epitopes? 
* Benchmarking
    * Local benchmarking tool to run our algorithm against others in FRED2
    * Add mhcflurry to FRED2
* Cancer-associated peptide databases
    * Look for association between peptides and cancer types
    * Add cancer peptides to bacteria/virus peptides and perform unsupervised learning
        * Since cancer peptides are just mutated human peptides, won't the clustering just separate based on species
        * Some peptides are only observed during development - may need to look at a developmental timecourse as a null set
* Viral sequence data set from Florian's collaborator
* Random peptide sequence generator (using RNN?)
    * Benchmark data set for MSA: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-66
    * Review of peptide simulation tools: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-184
