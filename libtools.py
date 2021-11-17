#!/usr/bin/env python3

import collections
from dict import *

def load(input):
    nucleotides = ['A', 'C', 'G', 'T']
    file = open(input[0])
    data = str(file.readlines()[input[1][0]:input[1][1]])
    genome = [x for x in data.upper() if (x in nucleotides)]
    genome = ''.join(genome)
    return genome

def translate(seq, dict):
    polypeptide = ""
    for i in range(0, len(seq)-3, 3):
        codon = seq[i:i + 3].replace('T', 'U')                                  # DNA to RNA transcription - Thymine is replaced with Uracil.
        amino = [k for k, v in dict.items() if codon in v]
        char = str(amino).split('/')[1].replace("']", "").strip()
        polypeptide += str(char)
    return polypeptide

def lookup_acid(acid):
    terminus = [k for k, v in mRNA_codon().items() if acid in k.split('/')[1]][0]
    return terminus

def amino_count(peptide):
    count = dict(collections.Counter(peptide))
    return count

def lookup_weight(peptide):
    water_mass = 18.01524
    weight = 0
    for i in peptide:
        weight += [v for k, v in molecular_weight().items() if i in k.split('/')[1]][0][0]
    weight -= water_mass * (len(peptide)-1)
    return weight

def lookup_halflife(acid):
    period = [v for k, v in halflife().items() if acid in k.split('/')[1]][0]
    return period

def hydropathy_index(peptide):
    index = 0
    for i in peptide:
         index += float([v for k, v in hydropathy().items() if i in k.split('/')[1]][0][0])
    index /= len(peptide)
    return index
