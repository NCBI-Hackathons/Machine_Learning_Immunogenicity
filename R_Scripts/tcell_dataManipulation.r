# Load required libraries
if (!require(plyr))
  install.packages("plyr")
if (!require(Peptides))
  install.packages("Peptides")

# Set working directory and read data
setwd("C:\\Pragyan\\Machine_Learning_Immunogenicity\\tcell_full_v3")

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


# aacomp(tcell[3,]$Epitope.Description)
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
tcell <- read.csv("tcell_full_v3.csv", header=F)
#tcell_orig <- read.csv("tcell_full_v3.csv", header=F) # Only for unit-testing; comment later
#tcell <- tcell_orig[1:25,] # Only for unit-testing; comment later

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

# Remove duplicates
tcell <- recode(tcell)
rownames(tcell) <- NULL

#Encode the peptides for Tensorflow using function peptides2input
tcell$Epitope.Description.calcValue <- sapply(1:nrow(tcell), function(x) peptides2input(as.character(tcell[x,]$Epitope.Description)))


#Use Peptide R package methods to calculate Indices and Theoretical Properties of Protein Sequences
# aacomp
calcValue <- data.frame(t(sapply(1:nrow(tcell), function(x) populatePeptide.aacomp(tcell[x,]$Epitope.Description))))
c1 <- sapply(1:ncol(calcValue), function(x) unlist(calcValue[,x]))
colnames(c1) <- colnames(calcValue)
tcell <- data.frame(cbind(tcell,c1))
rm(c1)

# aindex
Peptide.Aliphatic.Index <- sapply(1:nrow(tcell), function(x) aindex(tcell[x,]$Epitope.Description))
tcell <- cbind(tcell, Peptide.Aliphatic.Index)
# boman -- potential protein interaction
Peptide.Boman <- sapply(1:nrow(tcell), function(x) boman(tcell[x,]$Epitope.Description))
tcell <- cbind(tcell, Peptide.Boman)
# charge -- theoretical net charge of a protein sequence
Peptide.Charge <- sapply(1:nrow(tcell), function(x) charge(tcell[x,]$Epitope.Description))
tcell <- cbind(tcell, Peptide.Charge)
# hmoment
Peptide.hmoment <- sapply(1:nrow(tcell), function(x) hmoment(tcell[x,]$Epitope.Description))
tcell <- cbind(tcell, Peptide.hmoment)
# hydrophobicity
Peptide.hydrophobicity <- sapply(1:nrow(tcell), function(x) hydrophobicity(tcell[x,]$Epitope.Description))
tcell <- cbind(tcell, Peptide.hydrophobicity)
# instaindex
Peptide.instaindex <- sapply(1:nrow(tcell), function(x) instaindex(tcell[x,]$Epitope.Description))
tcell <- cbind(tcell, Peptide.instaindex)

# kidera with factor "helix.bend.pref"
Peptide.Kidera.helix.bend.pref <- sapply(1:nrow(tcell), function(x) kidera(tcell[x,]$Epitope.Description, "helix.bend.pref"))
tcell <- cbind(tcell, Peptide.Kidera.helix.bend.pref)
# kidera with factor "side.chain.size"
Peptide.Kidera.side.chain.size <- sapply(1:nrow(tcell), function(x) kidera(tcell[x,]$Epitope.Description, "side.chain.size"))
tcell <- cbind(tcell, Peptide.Kidera.side.chain.size)
# kidera with factor "extended.str.pref"
Peptide.Kidera.extended.str.pref <- sapply(1:nrow(tcell), function(x) kidera(tcell[x,]$Epitope.Description, "extended.str.pref"))
tcell <- cbind(tcell, Peptide.Kidera.extended.str.pref)
# kidera with factor "hydrophobicity"
Peptide.Kidera.hydrophobicity <- sapply(1:nrow(tcell), function(x) kidera(tcell[x,]$Epitope.Description, "hydrophobicity"))
tcell <- cbind(tcell, Peptide.Kidera.hydrophobicity)
# kidera with factor "double.bend.pref"
Peptide.Kidera.double.bend.pref <- sapply(1:nrow(tcell), function(x) kidera(tcell[x,]$Epitope.Description, "double.bend.pref"))
tcell <- cbind(tcell, Peptide.Kidera.double.bend.pref)
# kidera with factor "partial.spec.vol"
Peptide.Kidera.partial.spec.vol <- sapply(1:nrow(tcell), function(x) kidera(tcell[x,]$Epitope.Description, "partial.spec.vol"))
tcell <- cbind(tcell, Peptide.Kidera.partial.spec.vol)
# kidera with factor "flat.ext.pref"
Peptide.Kidera.flat.ext.pref <- sapply(1:nrow(tcell), function(x) kidera(tcell[x,]$Epitope.Description, "flat.ext.pref"))
tcell <- cbind(tcell, Peptide.Kidera.flat.ext.pref)
# kidera with factor "occurrence.alpha.reg"
Peptide.Kidera.occurrence.alpha.reg <- sapply(1:nrow(tcell), function(x) kidera(tcell[x,]$Epitope.Description, "occurrence.alpha.reg"))
tcell <- cbind(tcell, Peptide.Kidera.occurrence.alpha.reg)
# kidera with factor "pK.C"
Peptide.Kidera.pK.C <- sapply(1:nrow(tcell), function(x) kidera(tcell[x,]$Epitope.Description, "pK.C"))
tcell <- cbind(tcell, Peptide.Kidera.pK.C)
# kidera with factor "surrounding.hydropa"
Peptide.Kidera.surrounding.hydrop <- sapply(1:nrow(tcell), function(x) kidera(tcell[x,]$Epitope.Description, "surrounding.hydrop"))
tcell <- cbind(tcell, Peptide.Kidera.surrounding.hydrop)

# Ignore membpos
# mw
Peptide.mw <- sapply(1:nrow(tcell), function(x) mw(tcell[x,]$Epitope.Description))
tcell <- cbind(tcell, Peptide.mw)
# pI
Peptide.pI <- sapply(1:nrow(tcell), function(x) pI(tcell[x,]$Epitope.Description))
tcell <- cbind(tcell, Peptide.pI)


#Write out data for use in Tensorflow
write.csv(tcell,file="tcell_modified.csv",row.names=F)

#Calculate count by Organism
# count.Epitope <- count(tcell, vars=c("Epitope.Description", "Epitope.Organism.ID") )
# write.csv(count.Epitope,file="tcell_countByEpitopeOrganism.csv",row.names=F)


# Data value validation steps using histogram
#hist(nchar(as.character(bcell$Epitope.Description)))
#hist(log10(nchar(as.character(bcell$Epitope.Description))))
#quantile(nchar(as.character(bcell$Epitope.Description)))
#0%  25%  50%  75% 100% 
#2   15   15   15  829 
#quantile(nchar(as.character(bcell$Epitope.Description)),probs=1:10/10)
#10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
#  9   13   15   15   15   15   15   15   18  829 
#quantile(nchar(as.character(bcell$Epitope.Description)),probs=1:20/20)
#5%  10%  15%  20%  25%  30%  35%  40%  45%  50%  55%  60%  65%  70%  75%  80%  85%  90%  95% 100% 
#  8    9   10   13   15   15   15   15   15   15   15   15   15   15   15   15   15   18   21  829
# Pick Peptide length limit as 20