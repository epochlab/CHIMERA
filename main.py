#!/usr/bin/env python3

import zlib
import numpy as np
from libtools import *
from dict import mRNA_codon

SARS_CoV_2 = ('genome/SARS_CoV_2.txt', (692, 1191))
SARS_Tor2 = ('genome/SARS_Tor2.txt', (640, 1136))
MERS = ('genome/MERS.txt', (611, 1113))
HIV = ('genome/HIV.txt', (454, 608))
Ebola = ('genome/Ebola.txt', (371, 686))

genome = load(SARS_CoV_2)
# print(genome)

codon = mRNA_codon()

#Reading frame = START (ATG)
rf = genome.find('ATG')
if rf % 3==1:
    print('Reading frame:', rf)

print('Nucleobases:', len(genome))
print('C-G Content:', round((genome.count('C') + genome.count('G')) / len(genome)*100, 3), "%")
print('Compression:', len(zlib.compress(genome.encode("utf-8"))))

residue = translate(genome, codon)
print(residue.split('*'))

# >>> Identified proteins and amino acid count (SARS_CoV_2 only)
ORF1a = translate(genome[266-1: 13483], codon)                                         # ORF1a polyprotein - 4405
ORF1b = translate(genome[13468-1: 21555], codon)                                       # ORF1b polyprotein - 2695 overlapping sequence w/ ORF1a
S = translate(genome[21563-1: 25384], codon)                                           # Spike glycoprotein (structural) - 1273
ORF3a = translate(genome[25393-1: 26220], codon)                                       # ORF3a protein - 275
E = translate(genome[26245-1: 26472], codon)                                           # ORF4 envelope protein (structural) - 75
M = translate(genome[26523-1: 27191], codon)                                           # ORF5 membrane glycoprotein (structural) - 222
ORF6 = translate(genome[27202-1: 27387], codon)                                        # ORF6 protein - 61
ORF7a = translate(genome[27394-1: 27759], codon)                                       # ORF7a protein - 121
ORF7b = translate(genome[27756-1: 27887], codon)                                       # ORF7b protein - 43
ORF8 = translate(genome[27894-1: 28259], codon)                                        # ORF8 protein - 121
N = translate(genome[28274-1: 29533], codon)                                           # ORF9 nucleocapsid phosphoprotein (structural) - 419
ORF10 = translate(genome[29558-1: 29674], codon)                                       # ORF10 protein - 38
# print(ORF6)

index = 0
for pid, peptide in enumerate(residue.split('*')):
    if pid==index:

        n_terminus = lookup_acid(peptide[0])
        c_terminus = lookup_acid(peptide[-1])

        length = len(peptide)

        if length >= 2 and length <= 20:
            type = "Oligopeptide"
        else:
            type = "Polypeptide"

        mw = lookup_weight(peptide)
        hp = hydropathy_index(peptide)
        decay = lookup_halflife(peptide[0])
        aa_count = amino_count(peptide)
        formula, nb_atoms = atomic_composition(peptide)
        pos, neg = charged_residues(peptide)

        print("Chain Search:", residue.find(peptide))
        print(peptide)

        print("N-Terminus:", n_terminus, "| C-Terminus:", c_terminus)
        print("Sequence:", pid, "| Length:", length,  "| Type:", type, "| Molecular Weight (Da):", round(mw, 2), "| Half-life (N-end):", decay)
        print("Atomic Formula:", formula, "| Number of Atoms:", nb_atoms)
        print("Hydropathicity Index (GRAND Average):", round(hp, 3))

        print(aa_count)
        for a, c in aa_count.items():
            print(a, round(c * (100.0/length), 1), "%")

        print("+ charged residues (Arg | Lys | His):", pos)
        print("- charged residues (Asp | Glu):", neg)
