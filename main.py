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
    'End / STOP': ['UAA', 'UGA', 'UAG'],
    }

print("Amino acid keys:", len(dict.items()))

def translate(seq):
    polypeptide = ""
    for i in range(0, len(seq)-3, 3):
        codon = seq[i:i + 3].replace('T', 'U')
        amino = [k for k, v in dict.items() if codon in v]
        char = str(amino).split('/')[1].replace("']", "").strip()
        polypeptide += str(char)
    return polypeptide

pc = translate(genome)
print(pc.split('STOP'))

ORF1ab_polyprotein = (266, 13483)
spike_glycoprotein = (21563, 25384)
ORF6_protein = (27202, 27387)

gene_id = ORF6_protein

protein = translate(genome[gene_id[0]-1: gene_id[1]])
print(protein, len(protein))

index = "MFHLVDFQVTIAEILLIIMRTFKVSIWNLDYIINLIIKNLSKSLTENKYSQLDEEQPMEID"
print("Search:", index, pc.count(index))
