#!/usr/bin/env python3

nucleotides = ['A', 'C', 'G', 'T']

# Amino acid | Molecular weight | GRAVY | Codon
def mRNA_codon():
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
    return dict

def molecular_weight():
    dict = {
        'Ala / A': [89.09],
        'Ile / I': [131.18],
        'Arg / R': [174.20],
        'Leu / L': [131.18],
        'Asn / N': [132.12],
        'Lys / K': [146.19],
        'Asp / D': [133.10],
        'Met / M': [149.21],
        'Phe / F': [165.19],
        'Cys / C': [121.16],
        'Pro / P': [115.13],
        'Gln / Q': [146.15],
        'Ser / S': [105.09],
        'Glu / E': [147.13],
        'Thr / T': [119.12],
        'Trp / W': [204.23],
        'Gly / G': [75.07],
        'Tyr / Y': [181.19],
        'His / H': [155.16],
        'Val / V': [117.15],
        'STOP / *': [],
        }
    return dict

# HostID: Mammalian | Yeast | E.Coli
# Unit: Minutes
def halflife():
    dict = {
        'Ala / A': [264, 1200, 600],
        'Ile / I': [1200, 30, 600],
        'Arg / R': [60, 2, 2],
        'Leu / L': [330, 3, 2],
        'Asn / N': [84, 3, 600],
        'Lys / K': [78, 3, 2],
        'Asp / D': [66, 3, 600],
        'Met / M': [1800, 1200, 600],
        'Phe / F': [66, 3, 2],
        'Cys / C': [72, 1200, 600],
        'Pro / P': [1200, 1200, ""],
        'Gln / Q': [48, 10, 600],
        'Ser / S': [114, 1200, 600],
        'Glu / E': [60, 30, 600],
        'Thr / T': [432, 1200, 600],
        'Trp / W': [168, 3, 2],
        'Gly / G': [1800, 1200, 600],
        'Tyr / Y': [168, 10, 2],
        'His / H': [210, 10, 600],
        'Val / V': [6000, 1200, 600],
        'STOP / *': [],
        }
    return dict
