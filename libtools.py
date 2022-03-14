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

def reading_frame(seq):
    for i in range(3, len(seq)+1, 3):
        codon = seq[i-3:i]
        if codon == 'ATG':
            return i//3

def translate(seq, dict):
    polypeptide = ""
    for i in range(0, len(seq)-3, 3):
        codon = seq[i:i + 3].replace('T', 'U')                                  # DNA to RNA transcription - Thymine is replaced with Uracil.
        amino = [k for k, v in dict.items() if codon in v]
        char = str(amino).split('/')[1].replace("']", "").strip()
        polypeptide += str(char)
    return polypeptide

def lookup_value(input, dict):
    result = [v for k, v in dict.items() if input in k.split('/')[1]][0]
    return result

def lookup_amino(peptide):
    terminus = [k for k, v in mRNA_codon().items() if peptide in k.split('/')[1]][0]
    return terminus

def lookup_weight(peptide):
    water_mass = 18.01524
    weight = 0
    for i in peptide:
        weight += lookup_value(i, molecular_weight())
    weight -= water_mass * (len(peptide)-1)
    return weight

def lookup_halflife(peptide):
    period = lookup_value(peptide, halflife())
    return period

def hydropathy_index(peptide):
    index = 0
    for i in peptide:
         index += float(lookup_value(i, hydropathy()))
    index /= len(peptide)
    return index

def atomic_composition(peptide):
    chain = np.zeros((1,5), dtype=int)[0]
    for i in peptide:
        val = lookup_value(i, atomic())
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
    pos, neg = 0, 0
    for i in peptide:
        if i == "R" or i == "K" or i == "H":
            pos += 1
        if i == "D" or i == "E":
            neg += 1

    return pos, neg

def extinction_coefficient(peptide):
    nY, nW, nC = 0, 0, 0

    for i in peptide:
        if i == "Y":
            nY += 1
        if i == "W":
            nW += 1
        if i == "C":
            nC += 1

    # Ext. coefficient Tyrosine = 1490 | Tryptophan = 5500 | Cystine = 125
    # Cysteine does not absorb appreciably at wavelengths >260 nm, while Cystine does
    ext_coeff = (nY * 1490) + (nW * 5500) + (nC * 125)
    return ext_coeff

def instability_index(peptide):
    list = []
    for pid, amino in enumerate(peptide):
        layer = [v for k, v in DIWV().items() if amino in k][0]

        x = [v for k, v in layer.items() if amino in k][0]
        x_label = amino + amino

        if pid != len(peptide)-1:
            y = [v for k, v in layer.items() if peptide[pid+1] in k][0]
            y_label = amino + peptide[pid+1]
        else:
            y = 1.0
            y_label = 'NA'

        v = x * y
        list.append(v)

    II = (10/len(peptide)) * sum(list)
    return II

def aliphatic_index(peptide):
    nA, nV, nI, nL = 0, 0, 0, 0

    for i in peptide:
        if i == "A":
            nA += 1
        if i == "V":
            nV += 1
        if i == "I":
            nI += 1
        if i == "L":
            nL += 1

    # Aliphatic index = X(Ala) + a * X(Val) + b * ( X(Ile) + X(Leu) )
    # Coefficients A and B are the relative volume of the side chains (A = 2.9 | B = 3.9)
    total_atoms = len(peptide)
    index = nA + (2.9 * nV/total_atoms) + (3.9 * (nI/total_atoms + nL/total_atoms)) * 100
    return index
