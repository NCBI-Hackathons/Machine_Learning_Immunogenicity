# Load required libraries
if (!require(plyr))
  install.packages("plyr")
if (!require(Peptides))
  install.packages("Peptides")

# Set working directory and read data
setwd("C:\\Pragyan\\Machine_Learning_Immunogenicity\\tcell_full_v3")
tcell <- read.csv("tcell_full_v3.csv", header=F)

# Formulate header names
test <- tcell[1:2,]
names(tcell) <- apply(test,2,paste,collapse=" ")
names(tcell)
rm(test)
tcell <- tcell[-1,]
tcell <- tcell[-1,]

nrow(tcell) # 285725
ncol(tcell) # 95
names(tcell) <- make.names(names(tcell))

# Remove peptide types other than "Linear peptide" (field Epitope.Object.Type)
tcell <- tcell[tcell$Epitope.Object.Type == "Linear peptide",]

# Replace Positive-* in "Assay.Qualitative.Measure" with Positive
count(tcell, "Assay.Qualitative.Measure")
# Assay.Qualitative.Measure   freq
# 1                  Negative 168271
# 2                  Positive 102100
# 3             Positive-High   2433
# 4     Positive-Intermediate    937
# 5              Positive-Low   6992
tcell$Assay.Qualitative.Measure[tcell$Assay.Qualitative.Measure=="Positive-High"] = "Positive"
tcell$Assay.Qualitative.Measure[tcell$Assay.Qualitative.Measure=="Positive-Intermediate"] = "Positive"
tcell$Assay.Qualitative.Measure[tcell$Assay.Qualitative.Measure=="Positive-Low"] = "Positive"
count(tcell, "Assay.Qualitative.Measure")
# Assay.Qualitative.Measure   freq
# 1                  Negative 168271
# 2                  Positive 112462

#Remove non-peptides
tcell<- tcell[grepl("^[A-Z]+$",tcell$Epitope.Description),]
write.csv(tcell,file="tcell_modified.csv",row.names=F)

#Calculate count by Organism
count.Epitope <- count(tcell, vars=c("Epitope.Description", "Epitope.Organism.ID") )
write.csv(count.Epitope,file="tcell_countByEpitopeOrganism.csv",row.names=F)
