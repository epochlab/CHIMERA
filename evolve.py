#!/usr/bin/env python3

from libtools import load, translate
from codon import rna_table

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


pc = translate(genome, dict)
print(pc.split('*'))
