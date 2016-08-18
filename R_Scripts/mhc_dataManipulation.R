# Load required libraries
if (!require(plyr))
  install.packages("plyr")
if (!require(Peptides))
  install.packages("Peptides")
  
# Set working directory and read data
setwd("C:\\Pragyan\\Machine_Learning_Immunogenicity\\mhc_ligand_full")


# Functions
# For conversion of peptides to ascii num for use in Tensorflow
AA_ENCODING = list(
    A=2^0,
    C=2^1,
    D=2^2,
    E=2^3,
    F=2^4,
    G=2^5,
    H=2^6,
    I=2^7,
    K=2^8,
    L=2^9,
    M=2^10,
    N=2^11,
    P=2^12,
    Q=2^13,
    R=2^14,
    S=2^15,
    T=2^16,
    V=2^17,
    W=2^18,
    Y=2^19
)

peptides2input <- function(peptides) {
    AAs <- strsplit(peptides, '')
    unlist(lapply(AAs, function(x) {
        paste(sapply(x, function(y) AA_ENCODING[y]), collapse=',')
    }))
}


# aacomp(mhc[3,]$Epitope.Description)
# Number Mole%
# Tiny           4    20
# Small         10    50
# Aliphatic      4    20
# Aromatic       4    20
# NonPolar      12    60
# Polar          8    40
# Charged        5    25
# Basic          2    10
# Acidic         3    15
populatePeptide.aacomp <- function(peptide.value) {
  aacomp.retVal <- aacomp(peptide.value)
  
  retVal <- c()
  retVal$Tiny.Molecular.Percent <- c(0)
  retVal$Small.Molecular.Percent <- c(0)
  retVal$Aliphatic.Molecular.Percent <- c(0)
  retVal$Aromatic.Molecular.Percent <- c(0)
  retVal$Polar.Molecular.Percent <- c(0)
  retVal$Charged.Molecular.Percent <- c(0)
  retVal$Basic.Molecular.Percent <- c(0)
  retVal$Acidic.Molecular.Percent <- c(0)
  
  retVal$Tiny.Molecular.Percent <- aacomp.retVal[1,2]
  retVal$Small.Molecular.Percent <- aacomp.retVal[2,2]
  retVal$Aliphatic.Molecular.Percent <- aacomp.retVal[3,2]
  retVal$Aromatic.Molecular.Percent <- aacomp.retVal[4,2]
  retVal$Polar.Molecular.Percent <- aacomp.retVal[6,2]
  retVal$Charged.Molecular.Percent <- aacomp.retVal[7,2]
  retVal$Basic.Molecular.Percent <- aacomp.retVal[8,2]
  retVal$Acidic.Molecular.Percent <- aacomp.retVal[9,2]
  return (data.frame(retVal))
}


# Function to remove duplicates with conflicting Assay information
recode <- function(dataset) {
  # View(dataset[,94:ncol(dataset)])
  # Determine duplicates
   dataset$bin <- ifelse(dataset$Assay.Qualitative.Measure == 'Negative', 0, 1)
   peptides <- tapply(dataset$bin, dataset$Epitope.Description, function(x) c(sum(x==1), length(x)))
   pep1<-do.call(rbind, peptides)
   y<-as.numeric(pep1[,1])/as.numeric(pep1[,2])
   names(y) <- rownames(pep1)
   y<-y[y==0 | y == 1]
   w<-nchar(names(y)) <= 20
   y<-y[w]

   # Return a data frame of epitopes and assay measure only
   # This step is thus, also serving the dual purpose of removing unwanted columns of the dataset
   pep2<-data.frame(cbind(names(y), ifelse(y==0,0,1)))
   colnames(pep2) <- c("Epitope.Description","Assay.Qualitative.Measure")
   return (pep2)
}

# Read data

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

# Remove duplicates
mhc <- recode(mhc)
rownames(mhc) <- NULL



#Encode the peptides for Tensorflow using function peptides2input
mhc$Epitope.Description.calcValue <- sapply(1:nrow(mhc), function(x) peptides2input(as.character(mhc[x,]$Epitope.Description)))




#Use Peptide R package methods to calculate Indices and Theoretical Properties of Protein Sequences
# aacomp
calcValue <- data.frame(t(sapply(1:nrow(mhc), function(x) populatePeptide.aacomp(mhc[x,]$Epitope.Description))))
c1 <- sapply(1:ncol(calcValue), function(x) unlist(calcValue[,x]))
colnames(c1) <- colnames(calcValue)
tcell <- data.frame(cbind(mhc,c1))
rm(c1)


# aindex
Peptide.Aliphatic.Index <- sapply(1:nrow(mhc), function(x) aindex(mhc[x,]$Epitope.Description))
mhc <- cbind(mhc, Peptide.Aliphatic.Index)
# boman -- potential protein interaction
Peptide.Boman <- sapply(1:nrow(mhc), function(x) boman(mhc[x,]$Epitope.Description))
mhc <- cbind(mhc, Peptide.Boman)
# charge -- theoretical net charge of a protein sequence
Peptide.Charge <- sapply(1:nrow(mhc), function(x) charge(mhc[x,]$Epitope.Description))
mhc <- cbind(mhc, Peptide.Charge)
# hmoment
Peptide.hmoment <- sapply(1:nrow(mhc), function(x) hmoment(mhc[x,]$Epitope.Description))
mhc <- cbind(mhc, Peptide.hmoment)
# hydrophobicity
Peptide.hydrophobicity <- sapply(1:nrow(mhc), function(x) hydrophobicity(mhc[x,]$Epitope.Description))
mhc <- cbind(mhc, Peptide.hydrophobicity)
# instaindex
Peptide.instaindex <- sapply(1:nrow(mhc), function(x) instaindex(mhc[x,]$Epitope.Description))
mhc <- cbind(mhc, Peptide.instaindex)

# kidera with factor "helix.bend.pref"
Peptide.Kidera.helix.bend.pref <- sapply(1:nrow(mhc), function(x) kidera(mhc[x,]$Epitope.Description, "helix.bend.pref"))
mhc <- cbind(mhc, Peptide.Kidera.helix.bend.pref)
# kidera with factor "side.chain.size"
Peptide.Kidera.side.chain.size <- sapply(1:nrow(mhc), function(x) kidera(mhc[x,]$Epitope.Description, "side.chain.size"))
mhc <- cbind(mhc, Peptide.Kidera.side.chain.size)
# kidera with factor "extended.str.pref"
Peptide.Kidera.extended.str.pref <- sapply(1:nrow(mhc), function(x) kidera(mhc[x,]$Epitope.Description, "extended.str.pref"))
mhc <- cbind(mhc, Peptide.Kidera.extended.str.pref)
# kidera with factor "hydrophobicity"
Peptide.Kidera.hydrophobicity <- sapply(1:nrow(mhc), function(x) kidera(mhc[x,]$Epitope.Description, "hydrophobicity"))
mhc <- cbind(mhc, Peptide.Kidera.hydrophobicity)
# kidera with factor "double.bend.pref"
Peptide.Kidera.double.bend.pref <- sapply(1:nrow(mhc), function(x) kidera(mhc[x,]$Epitope.Description, "double.bend.pref"))
mhc <- cbind(mhc, Peptide.Kidera.double.bend.pref)
# kidera with factor "partial.spec.vol"
Peptide.Kidera.partial.spec.vol <- sapply(1:nrow(mhc), function(x) kidera(mhc[x,]$Epitope.Description, "partial.spec.vol"))
mhc <- cbind(mhc, Peptide.Kidera.partial.spec.vol)
# kidera with factor "flat.ext.pref"
Peptide.Kidera.flat.ext.pref <- sapply(1:nrow(mhc), function(x) kidera(mhc[x,]$Epitope.Description, "flat.ext.pref"))
mhc <- cbind(mhc, Peptide.Kidera.flat.ext.pref)
# kidera with factor "occurrence.alpha.reg"
Peptide.Kidera.occurrence.alpha.reg <- sapply(1:nrow(mhc), function(x) kidera(mhc[x,]$Epitope.Description, "occurrence.alpha.reg"))
mhc <- cbind(mhc, Peptide.Kidera.occurrence.alpha.reg)
# kidera with factor "pK.C"
Peptide.Kidera.pK.C <- sapply(1:nrow(mhc), function(x) kidera(mhc[x,]$Epitope.Description, "pK.C"))
mhc <- cbind(mhc, Peptide.Kidera.pK.C)
# kidera with factor "surrounding.hydropa"
Peptide.Kidera.surrounding.hydrop <- sapply(1:nrow(mhc), function(x) kidera(mhc[x,]$Epitope.Description, "surrounding.hydrop"))
mhc <- cbind(mhc, Peptide.Kidera.surrounding.hydrop)

# Ignore membpos
# mw
Peptide.mw <- sapply(1:nrow(mhc), function(x) mw(mhc[x,]$Epitope.Description))
mhc <- cbind(mhc, Peptide.mw)
# pI
Peptide.pI <- sapply(1:nrow(mhc), function(x) pI(mhc[x,]$Epitope.Description))
mhc <- cbind(mhc, Peptide.pI)

#Write out data for use in Tensorflow
write.csv(mhc,file="mhc_modified.csv",row.names=F)

#Calculate count by Organism
# count.bEpitope <- count(mhc, vars=c("Epitope.Description", "Epitope.Organism.ID") )
# write.csv(count.bEpitope,file="mhc_countByEpitopeOrganism.csv",row.names=F)
# Data value validation steps using histogram
#hist(nchar(as.character(mhc$Epitope.Description)))
#hist(log10(nchar(as.character(mhc$Epitope.Description))))
#quantile(nchar(as.character(mhc$Epitope.Description)))
#0%  25%  50%  75% 100% 
#2   15   15   15  829 
#quantile(nchar(as.character(mhc$Epitope.Description)),probs=1:10/10)
#10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
#  9   13   15   15   15   15   15   15   18  829 
#quantile(nchar(as.character(mhc$Epitope.Description)),probs=1:20/20)
#5%  10%  15%  20%  25%  30%  35%  40%  45%  50%  55%  60%  65%  70%  75%  80%  85%  90%  95% 100% 
#  8    9   10   13   15   15   15   15   15   15   15   15   15   15   15   15   15   18   21  829
# Pick Peptide length limit as 20
