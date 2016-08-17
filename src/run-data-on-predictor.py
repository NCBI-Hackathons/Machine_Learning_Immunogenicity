#!/usr/bin/env python
# coding: utf-8

# Author: F P Breitwieser

from __future__ import print_function, division, absolute_import
import numpy as np
import argparse
import sys
from sys import exit, stdin

from Fred2.Core import Allele, Peptide, Protein,generate_peptides_from_proteins
from Fred2.IO import read_lines, read_fasta
from Fred2.EpitopePrediction import EpitopePredictorFactory
#from nose.tools import eq_

import pepdata
import pandas as pd
pd.set_option('display.height', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 200)
pd.set_option('display.width', 200)


parser = argparse.ArgumentParser(description='Call epitope predictors on data.')
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--predictor', type=str, help='Epitope predictors [see all with --predictor=list]', required=True)
requiredNamed.add_argument('--dataset', type=str, help='Immunogenic dataset [see all with --dataset=list]', required=True)
parser.add_argument('-n', type=int, help='Number of rows to take from dataset')
parser.add_argument('--allele', type=str, help='HLA Type', default=["HLA-A*01:01","HLA-A*02:01","HLA-B*15:01"])

args = parser.parse_args()

all_predictors = [ name for name,version in EpitopePredictorFactory.available_methods().iteritems()]

all_predictors.remove("netmhcstabpan")
all_predictors.remove("netmhc")

if args.predictor == 'list':
	print("Set one of the predictors with --predictor:")
	print(all_predictors)
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
	_imm, _non = pepdata.iedb.tcell.load_classes(nrows=args.n)
elif args.dataset == 'iedb.mhc':
	_imm, _non = pepdata.iedb.mhc.load_classes(nrows=args.n)
elif args.dataset == 'imma2':
	_imm, _non = pepdata.imma2.load_classes()
elif args.dataset == 'stdin':
	_imm = [ line for line in stdin ]
else:
	print("available datasets: iedb.tcell, ideb.mhc, imma2")
	exit()

imm = [Peptide(elem) for elem in _imm]
non = [Peptide(elem) for elem in _non]

#hla='HLA-A2|HLA-A\*02'
#df = ()
#
#if args.dataset == 'iedb.tcell':
#	df =  pepdata.iedb.tcell.load_dataframe(nrows=args.n)
#elif args.dataset == 'iedb.mhc':
#	df =  pepdata.iedb.mhc.load_dataframe(nrows=args.n)
#elif args.dataset == 'iedb.bcell':
#	df =  pepdata.iedb.bcell.load_dataframe(nrows=args.n)
#elif args.dataset == 'stdin':
#	df = [ line for line in stdin ]
#else:
#	print("available datasets: iedb.tcell, ideb.mhc, imma2")
#	exit()
#
#print(df.head())
#

def run_predictor(pred, dataset):
	predictor = EpitopePredictorFactory(pred)
	results = ()
	try:
		results = predictor.predict(dataset, alleles=[ Allele(a) for a in args.allele ])
		print(results)
		print(results.describe())
	except ValueError:
		pass
	
	return(len(results),len(dataset))
	
def print_res(pred, res, dtype):
	print("RESULTS for %20s / %s: got scores for %5d/%d positives and %5d/%d negatives." % ((pred,dtype,) + res))

if args.predictor == 'all':
	for pred in all_predictors:
		try:
			print (pred," immunogenic peptides")
			print_res(pred, run_predictor(pred, imm), "    immunogenic peptides")
			if len(non) > 0:		
				print (pred," non-immunogenic peptides")
				print_res(pred, run_predictor(pred, non), "non-immunogenic peptides")
		except:
			pass
else:
	print (pred," immunogenic peptides")
	print_res(args.predictor,run_predictor(args.predictor, imm))
	if len(non) > 0:		
		print (pred," non-immunogenic peptides")
		print_res(args.predictor,run_predictor(args.predictor, non))


#print(results_imm)


#tcell_df = iedb.tcell.load_dataframe(nrows=1000)
#mhc_df = iedb.mhc.load_dataframe(nrows=1000)

## immunogenic and non-immunogenic sets
#imm, non = imma2.load_classes()


#peptides = [Peptide("SYFPEITHI"),Peptide("FIASNGVKL"), Peptide("LLGATCMFV")]
#allele = Allele("HLA-A*02:01")
#predictor = EpitopePredictorFactory("Syfpeithi")
#results = predictor.predict(peptides2, alleles=alleles)
#results.head()



