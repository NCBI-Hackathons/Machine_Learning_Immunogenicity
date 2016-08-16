library(plyr)

setwd("C:\\Pragyan\\Machine_Learning_Immunogenicity\\mhc_ligand_full")
mhc <- read.csv("mhc_ligand_full.csv", header=F)

# Formulate header names
test <- mhc[1:2,]
names(mhc) <- apply(test,2,paste,collapse=" ")
names(mhc)
rm(test)
mhc <- mhc[-1,]
mhc <- mhc[-1,]

nrow(mhc) # 510786
ncol(mhc) # 67
names(mhc) <- make.names(names(mhc))

# Remove peptide types other than "Linear peptide" (field Epitope.Object.Type)
mhc <- mhc[mhc$Epitope.Object.Type == "Linear peptide",]

# Replace Positive-* in "Assay.Qualitative.Measure" with Positive
count(mhc, "Assay.Qualitative.Measure")
# Assay.Qualitative.Measure   freq
# 1                  Negative 138709
# 2                  Positive 240742
# 3             Positive-High  48680
# 4     Positive-Intermediate  43308
# 5              Positive-Low  38955
mhc$Assay.Qualitative.Measure[mhc$Assay.Qualitative.Measure=="Positive-High"] = "Positive"
mhc$Assay.Qualitative.Measure[mhc$Assay.Qualitative.Measure=="Positive-Intermediate"] = "Positive"
mhc$Assay.Qualitative.Measure[mhc$Assay.Qualitative.Measure=="Positive-Low"] = "Positive"
count(mhc, "Assay.Qualitative.Measure")
# Assay.Qualitative.Measure   freq
# 1                  Negative 138709
# 2                  Positive 371685

#Remove non-peptides
mhc<- mhc[grepl("^[A-Z]+$",mhc$Epitope.Description),]

write.csv(mhc,file="mhc_modified.csv",row.names=F)


count.bEpitope <- count(mhc, vars=c("Epitope.Description", "Epitope.Organism.ID") )
write.csv(count.bEpitope,file="mhc_countByEpitopeOrganism.csv",row.names=F)
