#!/usr/bin/env python3

# Severe acute respiratory syndrome coronavirus 2 (COVID-19)
# SARS-CoV-2 is the virus that causes COVID-19 (coronavirus disease 2019), a positive-sense single-stranded RNA virus, that is contagious in humans and responsible for the COVID-19 pandemic.
# Dataset: NCBI RefSeq SARS-CoV-2 genome sequence record: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512

import sys
sys.path.append('.')

import zlib
from libtools import *
from dict import mRNA_codon

label, genome = load('data/NC_045512.2.txt')
# print(label, genome)

codon_table = mRNA_codon()

print(label.upper())
print('Nucleobases:', len(genome))
print('[START] Frame:', reading_frame(genome))
print('GC-Content:', round((genome.count('C') + genome.count('G')) / len(genome)*100, 3), "%")
print('Compression (zlib):', len(zlib.compress(genome.encode("utf-8"))))

# Sequential CGG position - Does NOT align with a modulo of 3, check reading_frame
print('CpG Islands:', genome.find('CGGCGG'))

res = translate(genome, codon_table)
print(res.split('*'))

# Reverse-engineered proteins
ORF1a = translate(genome[266-1: 13483], codon_table)                                         # ORF1a polyprotein - 4405
ORF1b = translate(genome[13468-1: 21555], codon_table)                                       # ORF1b polyprotein - 2695 overlapping sequence w/ ORF1a
S = translate(genome[21563-1: 25384], codon_table)                                           # Spike glycoprotein (structural) - 1273
ORF3a = translate(genome[25393-1: 26220], codon_table)                                       # ORF3a protein - 275
E = translate(genome[26245-1: 26472], codon_table)                                           # ORF4 envelope protein (structural) - 75
M = translate(genome[26523-1: 27191], codon_table)                                           # ORF5 membrane glycoprotein (structural) - 222
ORF6 = translate(genome[27202-1: 27387], codon_table)                                        # ORF6 protein - 61
ORF7a = translate(genome[27394-1: 27759], codon_table)                                       # ORF7a protein - 121
ORF7b = translate(genome[27756-1: 27887], codon_table)                                       # ORF7b protein - 43
ORF8 = translate(genome[27894-1: 28259], codon_table)                                        # ORF8 protein - 121
N = translate(genome[28274-1: 29533], codon_table)                                           # ORF9 nucleocapsid phosphoprotein (structural) - 419
ORF10 = translate(genome[29558-1: 29674], codon_table)                                       # ORF10 protein - 38
# print(S)

# Identify FURIN cleavage site in spike protein
print('FURIN cleavage site (Spike):', S.find('PRRAR'))
