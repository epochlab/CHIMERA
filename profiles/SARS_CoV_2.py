#!/usr/bin/env python3

# Severe acute respiratory syndrome coronavirus 2 (COVID-19)
# SARS-CoV-2 is the virus that causes COVID-19 (coronavirus disease 2019), a positive-sense single-stranded RNA virus, that is contagious in humans and responsible for the COVID-19 pandemic.
# Dataset: NCBI RefSeq SARS-CoV-2 genome sequence record: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512

# Immunoglobulin-M (IgM) Antibodies = Recent exposure
# IgG Antibodies = Created during initial infection but last for months

# June Almeida (1963)
# Coronavirus' are divided into 4 groups based on genetic difference
# Alpha, Beta, Gamma & Delta

# Alpha and Beta = Common
# Gamma and Delta = Rare

# Within the beta subdivision there are 3 different species:
# Lineage B
# SARS-like or SARS-related (SARSr)
# Sarbecoviruses

# Both SARS-1 and SARS-CoV-2 are both Sarbeviruses 
# Coronaviruses known to infect humans: 229E, NL63, HKU1, OC43, MERS, SARS-1, SARS-Cov-2

# Synonymous vs Non-synonymous
# Does NOT alter the residue chain = Synonymous
# DOES alter the residue chain = Non-Synonymous
# Caused by redundancy in the nucleotide sequence as some codons make the same protein.

# SARS-CoV-2 affects both upper and lower respiratory system.
# When in lungs and lower respiratory system.
# Affecting alveoli resulting in opaque shadows on scans.

# Severe reactions can trigger 'cytokine storm'.
# Inflammation caused by overreacting immune system leading to irreversible tissue damage.
# Affects other organs inc. intestines, heart, eyes, blood, sperm and CNS.
# The infection is systemic, causing damage to kidney, liver and spleen.

# Mojiang (Bengping) Copper Mine: N23°10'36' E101°21'28'
# https://s3.documentcloud.org/documents/6981198/Analysis-of-Six-Patients-With-Unknown-Viruses.pdf - Lu Xi, Masters Thesis (2013)
# 1,885km Bengpinghe > WIV



# -----

import sys
sys.path.append('../CHIMERA')

import zlib
from chimera.fasta import FastaIO
from chimera.measure import Signal
from chimera.codon import RNA

FASTA = FastaIO('NC_045512.2')

print("\n" + FASTA.label.upper())
print('Base Pairs:', len(FASTA.genome))
print('CG-Content:', round((FASTA.genome.count('C') + FASTA.genome.count('G')) / len(FASTA.genome)*100, 3), "%")
print('Compression (zlib):', Signal().compress(FASTA.genome))

# Sequential CGG position - Does NOT align with a modulo of 3, check reading_frame - ???
print('CpG Islands:', FASTA.genome.find('CGGCGG'))

print("\n" + FASTA.res)

# 15 genes encode 29 proteins
# 16 non-structural proteins which transform host cell into virus factory.
# nsp12, RNA-dependant RNA-polymerase (RdRp) - Copy / Generator Function
# nsp3, nsp4 and nsp6 recognise the internal structure of the host cell
# nsp14 proof reads the duplicate virus for error-checking. 

# Identify FURIN cleavage site in spike protein
print('FURIN cleavage site (Spike):', S.find('PRRAR'))

# Reverse-engineered proteins
ORF1a = FASTA.translate(FASTA.genome[266-1: 13483])                                         # ORF1a polyprotein - 4405
ORF1b = FASTA.translate(FASTA.genome[13468-1: 21555])                                       # ORF1b polyprotein - 2695 overlapping sequence w/ ORF1a
S = FASTA.translate(FASTA.genome[21563-1: 25384])                                           # Spike glycoprotein (structural) - 1273
ORF3a = FASTA.translate(FASTA.genome[25393-1: 26220])                                       # ORF3a protein - 275
E = FASTA.translate(FASTA.genome[26245-1: 26472])                                           # ORF4 envelope protein (structural) - 75
M = FASTA.translate(FASTA.genome[26523-1: 27191])                                           # ORF5 membrane glycoprotein (structural) - 222
ORF6 = FASTA.translate(FASTA.genome[27202-1: 27387])                                        # ORF6 protein - 61
ORF7a = FASTA.translate(FASTA.genome[27394-1: 27759])                                       # ORF7a protein - 121
ORF7b = FASTA.translate(FASTA.genome[27756-1: 27887])                                       # ORF7b protein - 43
ORF8 = FASTA.translate(FASTA.genome[27894-1: 28259])                                        # ORF8 protein - 121
N = FASTA.translate(FASTA.genome[28274-1: 29533])                                           # ORF9 nucleocapsid phosphoprotein (structural) - 419
ORF10 = FASTA.translate(FASTA.genome[29558-1: 29674])                                       # ORF10 protein - 38
# print(S)

# Angiotensin Converting Enzyme-II (ACE2)
# Receptor-Binding-Domain (RBD)
# ACE2 modulates blood pressure
