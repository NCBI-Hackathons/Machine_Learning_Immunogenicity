#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function, division, absolute_import
import numpy as np
import argparse
from sys import exit

from Fred2.Core import Allele, Peptide, Protein,generate_peptides_from_proteins
from Fred2.IO import read_lines, read_fasta
from Fred2.EpitopePrediction import EpitopePredictorFactory
#from nose.tools import eq_

from pepdata import iedb

parser = argparse.ArgumentParser(description='Call epitope predictors on data.')
parser.add_argument('predictor', type=str, help='Epitope predictors [see all with --predictor=list]')
parser.add_argument('dataset', type=str, help='Immunogenic dataset [see all with --dataset=list]')
parser.add_argument('-n', type=int, help='Number of rows to take from dataset')

args = parser.parse_args()


if args.predictor == 'list':
	print("Set one of the predictors with --predictor:")
	print([ name for name,version in EpitopePredictorFactory.available_methods().iteritems()].sort())
	print ("""
Details from https://bioinformatics.oxfordjournals.org/content/suppl/2016/02/26/btw113.DC1/S1.pdf
 SYFPEITHI     T-cell epitope  (Rammensee, et al., 1999)
 BIMAS         MHC-I binding   (Parker, et al., 1994)
 SVMHC         MHC-I binding   (Dönnes and Elofsson, 2002)
 ARB           MHC-I binding   (Bui, et al., 2005)
 SMM           MHC-I binding   (Peters and Sette, 2005)
 SMMPMBEC      MHC-I binding   (Kim, et al., 2009)
 Epidemix      MHC-I binding   (Feldhahn, et al., 2009)
 Comblib       MHC-I binding   (Sidney, et al., 2008)
 PickPocket*   MHC-I binding   (Zhang, et al., 2009)
 NetMHC*       MHC-I binding   (Lundegaard, et al., 2008)
 NetMHCpan*    MHC-I binding   (Hoof, et al., 2009)
 HAMMER        MHC-II binding  (Sturniolo, et al., 1999)
 TEPITOPEpan   MHC-II binding  (Zhang, et al., 2012)
 NetMHCII*     MHC-II binding  (Nielsen, et al., 2007)
 NetMHCIIpan*  MHC-II binding  (Karosiene, et al., 2013)
 UniTope       T-cell epitope  (Toussaint, et al., 2011)
 NetCTLpan*    T-cell epitope  (Stranzl, et al., 2010)
 """)

_imm, _non = set(), set()
if args.dataset == 'iedb.tcell':
	_imm, _non = iedb.tcell.load_classes(nrows=args.n)
elif args.dataset == 'iedb.bcell':
	_imm, _non = iedb.bcell.load_classes(nrows=args.n)
elif args.dataset == 'iedb.mhc':
	_imm, _non = iedb.mhc.load_classes(nrows=args.n)
else:
	print("available datasets: iedb.tcell, ideb.mhc")
	exit()

imm = [Peptide(elem) for elem in _imm]
non = [Peptide(elem) for elem in _non]

predictor = EpitopePredictorFactory(args.predictor)
results = predictor.predict(imm)
print("POSITIVES:")
print(results.head())

print("NEGATIVES:")
results = predictor.predict(non)
print(results.head())



#tcell_df = iedb.tcell.load_dataframe(nrows=1000)
#mhc_df = iedb.mhc.load_dataframe(nrows=1000)

## immunogenic and non-immunogenic sets
#imm, non = imma2.load_classes()


#peptides = [Peptide("SYFPEITHI"),Peptide("FIASNGVKL"), Peptide("LLGATCMFV")]
#allele = Allele("HLA-A*02:01")
#predictor = EpitopePredictorFactory("Syfpeithi")
#results = predictor.predict(peptides2, alleles=alleles)
#results.head()



