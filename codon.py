#!/usr/bin/env python3

nucleotides = ['A', 'C', 'G', 'T']

def rna_codon():
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
