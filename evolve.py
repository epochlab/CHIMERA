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

def lookup_acid(polypeptide, acid):
    terminus = [k for k, v in codon.items() if acid in k.split('/')[1]][0]
    return terminus

def amino_count(polypeptide):
    count = dict(collections.Counter(polypeptide))
    return count

def score(aa):
    return

def mutate(acid, epitope):
    return

def splice(sequence, epitope, mutant):
    return

def evolve(sequence, epitope):
    return

codon = mRNA_codon()

genome = load(('genome/SARS_CoV_2.txt', (692, 1191)))
residue = translate(genome, codon).split('*')

index = 5
water_atomicweight = 18.0153

for pid, polypeptide in enumerate(residue):
    if pid==index:

        length = len(polypeptide)
        n_terminus = lookup_acid(polypeptide, polypeptide[0])
        c_terminus = lookup_acid(polypeptide, polypeptide[-1])

        mw = 0.0
        for acid in polypeptide:
            mw += lookup_weight(acid)
        mw -= water_atomicweight * (length - 1)

        decay_rate = lookup_halflife(polypeptide[0])

        print("Sequence:", pid, "| Length:", length, "| Molecular Weight (Da):", mw, "| Half-life (N-end):", decay_rate)
        print(polypeptide)
        print("N-Terminus:", n_terminus, "| C-Terminus:", c_terminus)
        print(amino_count(polypeptide))
