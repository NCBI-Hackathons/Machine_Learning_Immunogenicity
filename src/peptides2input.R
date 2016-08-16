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
recode <- function(celltype, N) {
    tab<-readLines(gzfile(paste0(celltype, '.txt.gz')))
    tab<-strsplit(tab,"\t")
    tab<-lapply(tab, "[", 1:N)
    tab<-do.call(rbind, tab)
    colnames(tab) <- tab[1,]
    tab <- tab[-1,]
    tab <- as.data.frame(tab)
    tab$bin <- ifelse(tab$`Qualitative Measure` == 'Negative', 0, 1)
    tcell<-tab
    save(tcell, file=paste0(celltype, '.RData'))
    peptides <- tapply(tcell$bin, tcell$Description, function(x) c(sum(x==1), length(x)))
    pep1<-do.call(rbind, peptides)
    y<-as.numeric(pep1[,1])/as.numeric(pep1[,2])
    names(y) <- rownames(pep1)
    y<-y[y==0 | y == 1]
    w<-nchar(names(y)) <= 20
    y<-y[w]
    pep2<-peptides2input(names(y))
    pep3<-paste(pep2, ifelse(y==0,0,1), sep=',')
    writeLines(pep3, paste0(celltype, '_training_data_nodups.txt'))
    list(tab=tab,peptides=peptides,pep1=pep1,pep2=pep2,pep3=pep3)
}
