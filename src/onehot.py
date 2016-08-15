from bitarray import bitarray

AA_ENCODING = dict(
    A=bitarray('10000000000000000000'),
    C=bitarray('01000000000000000000'),
    D=bitarray('00100000000000000000'),
    E=bitarray('00010000000000000000'),
    F=bitarray('00001000000000000000'),
    G=bitarray('00000100000000000000'),
    H=bitarray('00000010000000000000'),
    I=bitarray('00000001000000000000'),
    K=bitarray('00000000100000000000'),
    L=bitarray('00000000010000000000'),
    M=bitarray('00000000001000000000'),
    N=bitarray('00000000000100000000'),
    P=bitarray('00000000000010000000'),
    Q=bitarray('00000000000001000000'),
    R=bitarray('00000000000000100000'),
    S=bitarray('00000000000000010000'),
    T=bitarray('00000000000000001000'),
    V=bitarray('00000000000000000100'),
    W=bitarray('00000000000000000010'),
    Y=bitarray('00000000000000000001')
)

class OneHotEncoded:
    def __init__(self, seq=None):
        """Create a new one-hot encoded protein sequence from a string."""
        if seq is not None:
            seq = seq.upper()
            self.size = len(seq)
            self.bits = bitarray()
            self.bits.encode(AA_ENCODING, seq)
    
    def __len__(self):
        return self.size
    
    def __getitem__(self, key):
        """Returns a new OneHotEncoded with the specified subsequence."""
        bits = self._subseq(key.start, key.stop)
        new_ohe = OneHotEncoded()
        new_ohe.size = len(bits) // 20
        new_ohe.bits = bits
        return new_ohe
    
    def copy(self, reverse_complement=False):
        new_ohe = OneHotEncoded()
        new_ohe.size = self.size
        new_ohe.bits = self.bits.copy()
        return new_ohe
    
    def _subseq(self, start, stop):
        return (
            self.bits[(start*20):(stop*20)],
            self.ambig[start:stop]
        )
    
    def to_bit_vector(self):
        return self.bits.to01()
            
    def to_aa(self):
        """Returns the sequence as a list of DNA base characters."""
        return self.bits.decode(BASE_ENCODING)
    
    def translate(self):
        """
        Returns a list of tuples where, at position i, the first element is the
        four-bit encoding of the AA at i, and the second element is the character
        corresponding to that base.
        """
        binstr = self.to_bit_vector()
        basestr = self.to_acgt()
        return ((binstr[(i*20):((i+1)*20)], basestr[i])
            for i in range(len(self)))
    
    def __repr__(self):
        """Returns the original DNA base string."""
        return "".join(self.to_acgt())
    
    def __str__(self):
        """
        Returns a string representation of the one-hot encoding in matrix form,
        along with the equivalent base character for each row.
        """
        return "ACGT\n----\n" + "\n".join("{} : {}".format(*t) for t in self.translate())
