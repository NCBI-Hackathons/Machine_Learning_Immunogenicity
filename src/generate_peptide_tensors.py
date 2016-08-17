from argparse import ArgumentParser
import numpy as np
import pandas as pd
from pickle import dump
from tqdm import tqdm

AA = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
AA_ENCODING = dict((a, 2**i) for i,a in enumerate(AA))

def generate_peptide_tensor(data, aaprops, output, maxlen = 20):
    peptides = pd.read_csv(data, sep="\t")
    npep = peptides.shape[0]
    aaprops = pd.read_csv(aaprops, sep="\t", index_col=2)
    aaprops = aaprops.iloc[:, 3:aaprops.shape[1]]
    nprop = aaprops.shape[1]
    shape = (npep, maxlen, nprop+1)
    tensor = np.zeros(shape)
    for i, p in tqdm(enumerate(peptides.iloc[:, 0])):
        if len(p) > maxlen:
            continue
        for j, a in enumerate(p):
            try:
                tensor[i, j, 0]  = AA_ENCODING[a]
                tensor[i, j, 1:] = aaprops.loc[a, :]
            except:
                print("Invalid AA: {} in {}".format(a, p))
    print("Writing tensor to file")
    with open(output, 'wb') as o:
        dump(tensor, o)
    
def main():
    parser = ArgumentParser()
    parser.add_argument("-a", "--aaprops", help="Amino Acid property file")
    parser.add_argument("-d", "--data", help="IEDB data file")
    parser.add_argument("-o", "--output", help="Output pickle file for Numpy array")
    args = parser.parse_args()
    
    generate_peptide_tensor(args.data, args.aaprops, args.output)

if __name__ == "__main__":
    main()
