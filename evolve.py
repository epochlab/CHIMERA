#!/usr/bin/env python3

import collections

from libtools import load, translate
from dict import mRNA_codon, molecular_weight, halflife

def lookup_weight(acid):
    weight = [v for k, v in molecular_weight().items() if acid in k.split('/')[1]][0][0]
    return weight

def lookup_halflife(acid):
    period = [v for k, v in halflife().items() if acid in k.split('/')[1]][0]
    return period

def amino_count(sequence):
    n = dict(collections.Counter(sequence))
    return n

def score(aa):
    return

def mutate(acid, epitope):
    return

def splice(sequence, epitope, mutant):
    return

def evolve(sequence, epitope):
    return

genome = load(('genome/SARS_CoV_2.txt', (692, 1191)))
residue = translate(genome, mRNA_codon()).split('*')

index = 0
water_atomicweight = 18.0153

for pid, polypeptide in enumerate(residue):
    if pid==index:

        length = len(polypeptide)

        mw = 0.0
        for acid in polypeptide:
            mw += lookup_weight(acid)
        mw -= water_atomicweight * (length - 1)

        n_terminus = polypeptide[0]
        c_terminus = polypeptide[-1]

        decay_rate = lookup_halflife(n_terminus)

        print("Sequence:", pid, "| Length:", length, "| Molecular Weight (Da):", mw, "| Half-life (N-end rule):", decay_rate)
        print(polypeptide)
        print(amino_count(polypeptide))
