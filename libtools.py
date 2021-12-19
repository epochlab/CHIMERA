#!/usr/bin/env python3

import numpy as np
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

def atomic_composition(peptide):
    chain = np.zeros((1,5), dtype=int)[0]
    for i in peptide:
        val = [v for k, v in atom().items() if i in k.split('/')[1]][0]
        chain = np.add(chain, val)

    atoms = ["C", "H", "N", "O", "S"]

    formula = []
    for tid, t in enumerate(chain):
        if tid == 1: # Hydrogen
            t -= (len(peptide)-1)*2
        if tid == 3: # Oxygen
            t -= (len(peptide)-1)
        group = atoms[tid]+str(t)
        formula.append(group)

    formula = ''.join(formula)
    nb_atoms = sum(chain) - (len(peptide)-1)*3

    return formula, nb_atoms

def amino_count(peptide):
    count = dict(collections.Counter(peptide))
    return count

def charged_residues(peptide):
    pos = 0
    neg = 0

    for i in peptide:
        if i == 'R' or i == 'K' or i == 'H':
            pos += 1
        if i == 'D' or i == 'E':
            neg += 1

    return pos, neg
