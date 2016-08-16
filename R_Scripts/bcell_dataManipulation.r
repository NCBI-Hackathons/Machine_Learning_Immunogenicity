library(plyr)

setwd("C:\\Pragyan\\Machine_Learning_Immunogenicity\\bcell_full_v3")
bcell <- read.csv("bcell_full_v3.csv", header=F)

# Formulate header names
test <- bcell[1:2,]
names(bcell) <- apply(test,2,paste,collapse=" ")
names(bcell)
rm(test)
bcell <- bcell[-1,]
bcell <- bcell[-1,]

nrow(bcell) # 385639
ncol(bcell) # 80
names(bcell) <- make.names(names(bcell))

# Remove peptide types other than "Linear peptide" (field Epitope.Object.Type)
bcell <- bcell[bcell$Epitope.Object.Type == "Linear peptide",]

# Replace Positive-* in "Assay.Qualitative.Measure" with Positive
count(bcell, "Assay.Qualitative.Measure")
# Assay.Qualitative.Measure   freq
# 1                  Negative 272246
# 2                  Positive  78545
# 3             Positive-High   2402
# 4     Positive-Intermediate    751
# 5              Positive-Low   5404
bcell$Assay.Qualitative.Measure[bcell$Assay.Qualitative.Measure=="Positive-High"] = "Positive"
bcell$Assay.Qualitative.Measure[bcell$Assay.Qualitative.Measure=="Positive-Intermediate"] = "Positive"
bcell$Assay.Qualitative.Measure[bcell$Assay.Qualitative.Measure=="Positive-Low"] = "Positive"
count(bcell, "Assay.Qualitative.Measure")
# Assay.Qualitative.Measure   freq
# 1                  Negative 272246
# 2                  Positive  87102

#Remove non-peptides
bcell<- bcell[grepl("^[A-Z]+$",bcell$Epitope.Description),]

write.csv(bcell,file="bcell_modified.csv",row.names=F)


count.bEpitope <- count(bcell, vars=c("Epitope.Description", "Epitope.Organism.ID") )
write.csv(count.bEpitope,file="bcell_countByEpitopeOrganism.csv",row.names=F)
