#!/usr/bin/env python3

nucleotides = ['A', 'C', 'G', 'T']

# Amino acid | Molecular weight | GRAVY | Codon
def rna_table():
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


# Molecular Weight | Half-life | Grand Average
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
