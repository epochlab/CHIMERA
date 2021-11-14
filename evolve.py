#!/usr/bin/env python3

from libtools import load, translate
from codon import rna_table, molecular_weight, halflife

def lookup_weight(acid):
    weight = [v for k, v in molecular_weight().items() if acid in k.split('/')[1]][0][0]
    return weight

def lookup_halflife(acid, hostid):
    weight = [v for k, v in halflife().items() if acid in k.split('/')[1]][0][hostid]
    return weight

def score(aa):
    return

def mutate(aa, epitope):
    return

def splice(sequence, epitope, mutant):
    return

def evolve(sequence, epitope):
    return

genome = load(('genome/SARS_CoV_2.txt', (692, 1191)))
residue = translate(genome, rna_table()).split('*')

index = 0

for pid, polypeptide in enumerate(residue):
    if pid==index:
        mw = 0.0
        hl = 0
        for acid in polypeptide:
            mw += lookup_weight(acid)

            val = lookup_halflife(acid, 0)
            if val>hl:
                hl=val

        print("Sequence:", pid, "| Length:", len(polypeptide), "| Molecular Weight:", mw, "| Half-life:", hl)
        print(polypeptide)
