#!/usr/bin/env python3

nucleotides = ['A', 'C', 'G', 'T']

# Amino acid | Molecular weight | Codon
def rna_table():
    dict = {
        'Ala / A': ['89.09', 'GCU', 'GCC', 'GCA', 'GCG'],
        'Ile / I': ['131.18', 'AUU', 'AUC', 'AUA'],
        'Arg / R': ['174.20', 'CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'Leu / L': ['131.18', 'CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG'],
        'Asn / N': ['132.12', 'AAU', 'AAC'],
        'Lys / K': ['146.19', 'AAA', 'AAG'],
        'Asp / D': ['133.10', 'GAU', 'GAC'],
        'Met / M': ['149.21', 'AUG'],
        'Phe / F': ['165.19', 'UUU', 'UUC'],
        'Cys / C': ['121.16', 'UGU', 'UGC'],
        'Pro / P': ['115.13', 'CCU', 'CCC', 'CCA', 'CCG'],
        'Gln / Q': ['146.15', 'CAA', 'CAG'],
        'Ser / S': ['105.09', 'UCU', 'UCC', 'UCA', 'UCG',  'AGU', 'AGC'],
        'Glu / E': ['147.13', 'GAA', 'GAG'],
        'Thr / T': ['119.12', 'ACU', 'ACC', 'ACA', 'ACG'],
        'Trp / W': ['204.23', 'UGG'],
        'Gly / G': ['75.07', 'GGU', 'GGC', 'GGA', 'GGG'],
        'Tyr / Y': ['181.19', 'UAU', 'UAC'],
        'His / H': ['155.16', 'CAU', 'CAC'],
        'Val / V': ['117.15', 'GUU', 'GUC', 'GUA', 'GUG'],
        'STOP / *': ['UAA', 'UGA', 'UAG'],
        }
    return dict
