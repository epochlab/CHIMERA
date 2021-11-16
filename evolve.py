#!/usr/bin/env python3

import collections

from libtools import load, translate
from dict import mRNA_codon, molecular_weight, halflife

def lookup_acid(acid):
    terminus = [k for k, v in codon.items() if acid in k.split('/')[1]][0]
    return terminus

def lookup_weight(acid):
    weight = [v for k, v in molecular_weight().items() if acid in k.split('/')[1]][0][0]
    return weight

def lookup_halflife(acid):
    period = [v for k, v in halflife().items() if acid in k.split('/')[1]][0]
    return period

def amino_count(peptide):
    count = dict(collections.Counter(peptide))
    return count

def score(aa):
    return

def mutate(acid, epitope):
    return

def splice(sequence, epitope, mutant):
    return

def evolve(sequence, epitope):
    return

water_mass = 18.01524

codon = mRNA_codon()

genome = load(('genome/SARS_CoV_2.txt', (692, 1191)))
residue = translate(genome, codon)

index = 5
for pid, peptide in enumerate(residue.split('*')):
    if pid==index:

        n_terminus = lookup_acid(peptide[0])
        c_terminus = lookup_acid(peptide[-1])

        length = len(peptide)

        if length >= 2 and length <= 20:
            type = "Oligopeptide"
        else:
            type = "Polypeptide"

        mw = 0.0
        for acid in peptide:
            mw += lookup_weight(acid)
        mw -= water_mass * (length-1)

        decay_rate = lookup_halflife(peptide[0])
        composition = amino_count(peptide)

        print(peptide)
        print("N-Terminus:", n_terminus, "| C-Terminus:", c_terminus)
        print("Sequence:", pid, "| Length:", length,  "| Type:", type, "| Molecular Weight (Da):", mw, "| Half-life (N-end):", decay_rate)
        print("Chain Search:", residue.find(peptide))
        print(composition)

        for a, c in composition.items():
            print(a + ':', "{:.3f}".format(c * (100.0/length)), "%")
