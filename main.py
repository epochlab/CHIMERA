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

mRNA_dict = {
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

print("Amino acid keys:", len(mRNA_dict.items()))

def decode(seq):
    protein = ""
    for i in range(0, len(seq)-3, 3):
        codon = seq[i:i + 3].replace('T', 'U')
        amino = [k for k, v in mRNA_dict.items() if codon in v]
        char = str(amino).split('/')[1].replace("']", "").strip()
        protein += str(char)
    return protein

prt = decode(genome)
print(prt.split('STOP'))

index = "MFHLVDFQVTIAEILLIIMRTFKVSIWNLDYIINLIIKNLSKSLTENKYSQLDEEQPMEID"
print("Search:", index, prt.count(index))

gene_id = (27202, 27387)
ORF1ab = decode(genome[gene_id[0]-1: gene_id[1]])
print(ORF1ab)
