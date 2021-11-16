#!/usr/bin/env python3

import zlib
from libtools import load, translate
from dict import mRNA_codon

SARS_CoV_2 = ('genome/SARS_CoV_2.txt', (692, 1191))
SARS_Tor2 = ('genome/SARS_Tor2.txt', (640, 1136))
MERS = ('genome/MERS.txt', (611, 1113))
HIV = ('genome/HIV.txt', (454, 608))
Ebola = ('genome/Ebola.txt', (371, 686))

genome = load(SARS_CoV_2)
# print(genome)

#DNA: ACGT | RNA: ACGU
codon = mRNA_codon()

#Reading frame = START (ATG)
rf = genome.find('ATG')
if rf % 3==1:
    print('Reading frame:', rf)

print('Nucleobases:', len(genome))
print('C-G Content:', "{:.3f}".format(((genome.count('C') + genome.count('G')) / len(genome)*100)), "%")
print('Compression:', len(zlib.compress(genome.encode("utf-8"))))

residue = translate(genome, codon)
# print(residue.split('*'))

### Identified proteins and amino acid count (SARS_CoV_2 only)
ORF1a = translate(genome[266-1: 13483], codon)                                         # ORF1a polyprotein - 4405
ORF1b = translate(genome[13468-1: 21555], codon)                                       # ORF1b polyprotein - 2695 overlapping sequence w/ ORF1a
S = translate(genome[21563-1: 25384], codon)                                           # Spike glycoprotein (structural) - 1273
ORF3a = translate(genome[25393-1: 26220], codon)                                       # ORF3a protein - 275
E = translate(genome[26245-1: 26472], codon)                                           # ORF4 envelope protein (structural) - 75
M = translate(genome[26523-1: 27191], codon)                                           # ORF5 membrane glycoprotein (structural) - 222
ORF6 = translate(genome[27202-1: 27387], codon)                                        # ORF6 protein - 61
ORF7a = translate(genome[27394-1: 27759], codon)                                       # ORF7a protein - 121
ORF7b = translate(genome[27756-1: 27887], codon)                                       # ORF7b protein - 43
ORF8 = translate(genome[27894-1: 28259], codon)                                        # ORF8 protein - 121
N = translate(genome[28274-1: 29533], codon)                                           # ORF9 nucleocapsid phosphoprotein (structural) - 419
ORF10 = translate(genome[29558-1: 29674], codon)                                       # ORF10 protein - 38

# print(ORF6)
