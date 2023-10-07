#!/usr/bin/env python3

# Immunoglobulin-M (IgM) Antibodies = Recent exposure
# IgG Antibodies = Created during initial infection but last for months

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

# -----

import sys
sys.path.append('..')

import zlib
from chimera.fasta import FastaIO
from chimera.measure import Molecule

FASTA = FastaIO('NC_045512.2')

# label, genome = FASTA.load('../genome/NC_045512.2.fasta')
# print(label, genome)

# print(FASTA.label.upper())
# print('Base Pairs:', len(M.genome))
# print('CG-Content:', round((M.genome.count('C') + M.genome.count('G')) / len(M.genome)*100, 3), "%")
# print('Compression (zlib):', len(zlib.compress(M.genome.encode("utf-8"))))

# # Sequential CGG position - Does NOT align with a modulo of 3, check reading_frame - ???
# print('CpG Islands:', M.genome.find('CGGCGG'))

# # print(M.res)

# # 15 genes encode 29 proteins
# # 16 non-structural proteins which transform host cell into virus factory.
# # nsp12, RNA-dependant RNA-polymerase (RdRp) - Copy / Generator Function
# # nsp3, nsp4 and nsp6 recognise the internal structure of the host cell
# # nsp14 proof reads the duplicate virus for error-checking. 

# # Reverse-engineered proteins
# ORF1a = M._translate(M.genome[266-1: 13483], RNA())                                         # ORF1a polyprotein - 4405
# ORF1b = M._translate(M.genome[13468-1: 21555], RNA())                                       # ORF1b polyprotein - 2695 overlapping sequence w/ ORF1a
# S = M._translate(M.genome[21563-1: 25384], RNA())                                           # Spike glycoprotein (structural) - 1273
# ORF3a = M._translate(M.genome[25393-1: 26220], RNA())                                       # ORF3a protein - 275
# E = M._translate(M.genome[26245-1: 26472], RNA())                                           # ORF4 envelope protein (structural) - 75
# M = M._translate(M.genome[26523-1: 27191], RNA())                                           # ORF5 membrane glycoprotein (structural) - 222
# ORF6 = M._translate(M.genome[27202-1: 27387], RNA())                                        # ORF6 protein - 61
# ORF7a = M._translate(M.genome[27394-1: 27759], RNA())                                       # ORF7a protein - 121
# ORF7b = M._translate(M.genome[27756-1: 27887], RNA())                                       # ORF7b protein - 43
# ORF8 = M._translate(M.genome[27894-1: 28259], RNA())                                        # ORF8 protein - 121
# N = M._translate(M.genome[28274-1: 29533], RNA())                                           # ORF9 nucleocapsid phosphoprotein (structural) - 419
# ORF10 = M._translate(M.genome[29558-1: 29674], RNA())                                       # ORF10 protein - 38
# print(S)

# # Angiotensin Converting Enzyme-II (ACE2)
# # Receptor-Binding-Domain (RBD)
# # ACE2 modulates blood pressure
