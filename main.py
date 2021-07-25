#!/usr/bin/env python3

import zlib

nucleotides = ['A', 'C', 'G', 'T']

file = open('genome.txt')

data = str(file.readlines()[692:1191])
genome = [x for x in data.upper() if (x in nucleotides)]
genome = ''.join(genome)

length = len(genome)
dst = zlib.compress(genome.encode("utf-8"))

print(genome)
print('Base Pairs:', length)
print('Compression:', len(dst))

dict = {
    'Ala / A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'Ile / I': ['AUU', 'AUC', 'AUA'],
    'Arg / R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'Leu / L': ['CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG'],
    'Asn / N': ['AAU', 'AAC'],
    'Lys / K': ['AAA', 'AAG'],
    'Asp / D': ['GAU', 'GAC'],
    'Met / M': ['AUG'],
    'Phe / F': ['UUU', 'UUC'],
    'Cys / C': ['UGU', 'UGC'],
    'Pro / P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Gln / Q': ['CAA', 'CAG'],
    'Ser / S': ['UCU', 'UCC', 'UCA', 'UCG',  'AGU', 'AGC'],
    'Glu / E': ['GAA', 'GAG'],
    'Thr / T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'Trp / W': ['UGG'],
    'Gly / G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'Tyr / Y': ['UAU', 'UAC'],
    'His / H': ['CAU', 'CAC'],
    'Val / V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'STOP / *': ['UAA', 'UGA', 'UAG'],
    }

print("Amino acids:", len(dict.items()))

def translate(seq):
    polypeptide = ""
    for i in range(0, len(seq)-3, 3):
        codon = seq[i:i + 3].replace('T', 'U')
        amino = [k for k, v in dict.items() if codon in v]
        char = str(amino).split('/')[1].replace("']", "").strip()
        polypeptide += str(char)
    return polypeptide

pc = translate(genome)
print(pc.split('*'))

index = "MFHLVDFQVTIAEILLIIMRTFKVSIWNLDYIINLIIKNLSKSLTENKYSQLDEEQPMEID"
print("Search:", index, pc.count(index))

ORF1a = translate(genome[266-1: 13483])                                         # ORF1a polyprotein - 4405
ORF1b = translate(genome[13468-1: 21555])                                       # ORF1b polyprotein - 2695*
S = translate(genome[21563-1: 25384])                                           # Spike glycoprotein (structural) - 1273
ORF3a = translate(genome[25393-1: 26220])                                       # ORF3a protein - 275
E = translate(genome[26245-1: 26472])                                           # ORF4 envelope protein (structural) - 75
M = translate(genome[26523-1: 27191])                                           # ORF5 membrane glycoprotein (structural) - 222
ORF6 = translate(genome[27202-1: 27387])                                        # ORF6 protein - 61
ORF7a = translate(genome[27394-1: 27759])                                       # ORF7a protein - 121
ORF7b = translate(genome[27756-1: 27887])                                       # ORF7b protein - 43
ORF8 = translate(genome[27894-1: 28259])                                        # ORF8 protein - 121
N = translate(genome[28274-1: 29533])                                           # ORF9 nucleocapsid phosphoprotein (structural) - 419
ORF10 = translate(genome[29558-1: 29674])                                       # ORF10 protein - 38

disp = ORF1a
print(disp, len(disp))
