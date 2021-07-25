#!/usr/bin/env python3

import zlib
from libtools import load, translate
from codon import rna_codon

SARS_CoV_2 = ('genome/SARS_CoV_2.txt', (692, 1191))
SARS_Tor2 = ('genome/SARS_Tor2.txt', (640, 1136))
MERS = ('genome/MERS.txt', (611, 1113))
HIV = ('genome/HIV.txt', (454, 608))
Ebola = ('genome/Ebola_sudan.txt', (371, 686))

genome = load(SARS_CoV_2)

# print(genome)
print('Base pairs:', len(genome))
print('G-C content:', ((genome.count('G') + genome.count('C')) / len(genome) * 100))

rf = genome.find('ATG')
if rf % 3==1:
    print('Reading frame:', rf)

print('Compression:', len(zlib.compress(genome.encode("utf-8"))))

dict = rna_codon()
print("Amino acids:", len(dict.items()))

pc = translate(genome, dict)
# print(pc.split('*'))

index = "MFHLVDFQVTIAEILLIIMRTFKVSIWNLDYIINLIIKNLSKSLTENKYSQLDEEQPMEID"
print("Chain search:", index, pc.find(index))

# ORF1a = translate(genome[266-1: 13483], dict)                                         # ORF1a polyprotein - 4405
# ORF1b = translate(genome[13468-1: 21555], dict)                                       # ORF1b polyprotein - 2695 overlap sequence
# S = translate(genome[21563-1: 25384], dict)                                           # Spike glycoprotein (structural) - 1273
# ORF3a = translate(genome[25393-1: 26220], dict)                                       # ORF3a protein - 275
# E = translate(genome[26245-1: 26472], dict)                                           # ORF4 envelope protein (structural) - 75
# M = translate(genome[26523-1: 27191], dict)                                           # ORF5 membrane glycoprotein (structural) - 222
# ORF6 = translate(genome[27202-1: 27387], dict)                                        # ORF6 protein - 61
# ORF7a = translate(genome[27394-1: 27759], dict)                                       # ORF7a protein - 121
# ORF7b = translate(genome[27756-1: 27887], dict)                                       # ORF7b protein - 43
# ORF8 = translate(genome[27894-1: 28259], dict)                                        # ORF8 protein - 121
# N = translate(genome[28274-1: 29533], dict)                                           # ORF9 nucleocapsid phosphoprotein (structural) - 419
# ORF10 = translate(genome[29558-1: 29674], dict)                                       # ORF10 protein - 38
#
# # disp = S
# # print(disp, len(disp))
