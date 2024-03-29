#!/usr/bin/env python3

import numpy as np
import collections, math, zlib
from PIL import Image
from chimera.codon import *

def lookup_value(input, dict):
    return [v for k, v in dict.items() if input in k.split('/')[1]][0]

def lookup_amino(peptide):
    return [k for k, _ in RNA().items() if peptide in k.split('/')[1]][0]

def gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq)*100

def peptype(peptide):
    length = len(peptide)
    if length >= 2 and length <= 20:
        type = "Oligopeptide"
    elif length > 20:
        type = "Polypeptide"
    return type

def molecular_weight(peptide):
    water_mass = 18.01524
    weight = 0
    for i in peptide:
        weight += lookup_value(i, amino_weight())
    weight -= water_mass * (len(peptide)-1)
    return weight

def charge_at_pH(pH, peptide):
    alpha_amino = [v for k, v in pKa().items() if k == 'alpha_amino'][0]
    alpha_carboxy = [v for k, v in pKa().items() if k == 'alpha_carboxy'][0]
    sidechain_positive = [v for k, v in pKa().items() if k == 'sidechain_positive'][0]
    sidechain_negative = [v for k, v in pKa().items() if k == 'sidechain_negative'][0]
    nterm, cterm = peptide[0], peptide[-1]

    net_charge = 0.0
    net_charge += math.pow(10, alpha_amino[nterm]) / (math.pow(10, alpha_amino[nterm]) + math.pow(10, pH))
    net_charge -= math.pow(10, pH) / (math.pow(10, alpha_carboxy[cterm]) + math.pow(10, pH))

    for aa in peptide:
        if aa in sidechain_positive:
            net_charge += math.pow(10, sidechain_positive[aa]) / (math.pow(10, sidechain_positive[aa]) + math.pow(10, pH))
        if aa in sidechain_negative:
            net_charge -= math.pow(10, pH) / (math.pow(10, sidechain_negative[aa]) + math.pow(10, pH))

    return net_charge

def isoelectric_point(peptide, pH=7.0, min=4, max=12):
    charge = charge_at_pH(pH, peptide)
    if max - min > 0.0001:
        if charge > 0.0:
            min = pH
        else:
            max = pH
        next_pH = (min + max) / 2
        return isoelectric_point(peptide, next_pH, min, max)
    return pH

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
    return dict(collections.Counter(peptide))

def charged_residues(peptide):
    pos, neg = 0, 0
    for i in peptide:
        if i == "R" or i == "K" or i == "H":
            pos += 1
        if i == "D" or i == "E" or i == "C" or i == "Y":
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
        if pid != len(peptide)-1:
            val = [v for k, v in layer.items() if peptide[pid+1] in k][0]
            label = amino + peptide[pid+1]
        else:
            val = 0.0
            label = 'NA'
        list.append(val)

    II = (10/(len(peptide))) * sum(list)
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

def compress(seq):
    return len(zlib.compress(seq.encode('utf-8')))

def binary_array_to_hex(arr):
    bit_string = ''.join(str(b) for b in 1 * arr.flatten())
    width = int(np.ceil(len(bit_string)/4))
    return '{:0>{width}x}'.format(int(bit_string, 2), width=width)

def seq_to_pixels(seq):
    pixels = []

    for n in seq:
        if n == 'A':
            pixels.append((30,30,30))
        if n == 'G':
            pixels.append((65,150,65))
        if n == 'C':
            pixels.append((255,255,255))
        if n == 'T':
            pixels.append((0,50,140))

    dim = round(math.sqrt(len(seq)))
    diff = dim ** 2 - len(seq)

    for _ in range(diff):
        pixels.append((255,255,0))

    array = np.array(pixels, dtype=np.uint8)
    return Image.frombytes("RGB", (dim, dim), bytes(array))

def average_hash(image, hash_size=8, mean=np.mean):
    img = image.convert('L').resize((hash_size, hash_size), Image.ANTIALIAS)
    pixels = np.asarray(img)
    return binary_array_to_hex(pixels > mean(pixels))
