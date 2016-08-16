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
