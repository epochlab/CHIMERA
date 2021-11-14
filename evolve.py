#!/usr/bin/env python3

from libtools import load, translate
from codon import rna_table

def molecularweight(acid, dict):
    lookup = float([v for k, v in dict.items() if acid in k.split('/')[1]][0][0])
    return lookup

def score(aa):
    return

def mutate(aa, epitope):
    return

def splice(sequence, epitope, mutant):
    return

def evolve(sequence, epitope):
    return

dict = rna_table()

genome = load(('genome/SARS_CoV_2.txt', (692, 1191)))
residue = translate(genome, dict).split('*')

index = 12

for pid, polypeptide in enumerate(residue):
    if pid==index:
        mw = 0.0
        for char in polypeptide:
            mw += molecularweight(char, dict)

        print("Sequence:", pid, "| Length:", len(polypeptide), "| Molecular Weight:", mw)
        print(polypeptide)
