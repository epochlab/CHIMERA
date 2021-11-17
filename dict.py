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
        'Pro / P': [1200, 1200, 600],
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

def hydropathy():
    dict = {
        'Ala / A': [1.8],
        'Ile / I': [4.5],
        'Arg / R': [-4.5],
        'Leu / L': [3.8],
        'Asn / N': [-3.5],
        'Lys / K': [-3.9],
        'Asp / D': [-3.5],
        'Met / M': [1.9],
        'Phe / F': [2.8],
        'Cys / C': [2.5],
        'Pro / P': [-1.6],
        'Gln / Q': [-3.5],
        'Ser / S': [-0.8],
        'Glu / E': [-3.5],
        'Thr / T': [-0.7],
        'Trp / W': [-0.9],
        'Gly / G': [-0.4],
        'Tyr / Y': [-1.3],
        'His / H': [-3.2],
        'Val / V': [4.2],
        'STOP / *': [],
        }
    return dict

def atom():
    dict = {
        'Ala / A': {C: 3, H: 7 , O: 2, N: 1, S: 0}, # C3H7NO2
        'Ile / I': {C: 6, H: 13, O: 2, N: 1, S: 0}, # C6H13NO2
        'Arg / R': {C: 6, H: 14, O: 2, N: 4, S: 0}, # C6H14N4O2
        'Leu / L': {C: 6, H: 13, O: 2, N: 1, S: 0}, # C6H13NO2
        'Asn / N': {C: 4, H: 8 , O: 3, N: 2, S: 0}, # C4H8N2O3
        'Lys / K': {C: 6, H: 14, O: 2, N: 2, S: 0}, # C6H14N2O2
        'Asp / D': {C: 4, H: 7 , O: 4, N: 1, S: 0}, # C4H7NO4
        'Met / M': {C: 5, H: 11 ,O: 2, N: 1, S: 1}, # C5H11NO2S
        'Phe / F': {C: 9, H: 11, O: 2, N: 1, S: 0}, # C9H11NO2
        'Cys / C': {C: 3, H: 7 , O: 2, N: 1, S: 1}, # C3H7NO2S
        'Pro / P': {C: 5, H: 9 , O: 2, N: 1, S: 0}, # C5H9NO2
        'Gln / Q': {C: 5, H: 10, O: 3, N: 2, S: 0}, # C5H10N2O3
        'Ser / S': {C: 3, H: 7 , O: 3, N: 1, S: 0}, # C3H7NO3
        'Glu / E': {C: 5, H: 9 , O: 4, N: 1, S: 0}, # C5H9NO4
        'Thr / T': {C: 4, H: 9 , O: 3, N: 1, S: 0}, # C4H9NO3
        'Trp / W': {C: 11,H: 12, O: 2, N: 2, S: 0}, # C11H12N2O2
        'Gly / G': {C: 2, H: 5 , O: 2, N: 1, S: 0}, # C2H5NO2
        'Tyr / Y': {C: 9, H: 11, O: 3, N: 1, S: 0}, # C9H11NO3
        'His / H': {C: 6, H: 9 , O: 2, N: 3, S: 0}, # C6H9N3O2
        'Val / V': {C: 5, H: 11, O: 2, N: 1, S: 0}, # C5H11NO2
        'STOP / *': [],
        }
    return dict
