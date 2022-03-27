#!/usr/bin/env python3

import zlib
from libtools import *
from dict import mRNA_codon

SARS_CoV_2 = ('genome/SARS_CoV_2.txt', (692, 1191))
SARS_Tor2 = ('genome/SARS_Tor2.txt', (640, 1136))
MERS = ('genome/MERS.txt', (611, 1113))
HIV = ('genome/HIV.txt', (454, 608))
Ebola = ('genome/Ebola.txt', (371, 686))

genome = load(SARS_CoV_2)
# print(genome)

codon_table = mRNA_codon()

print('Nucleobases:', len(genome))
print('[START] Frame:', reading_frame(genome))
print('GC-Content:', round((genome.count('C') + genome.count('G')) / len(genome)*100, 3), "%")
print('Compression (zlib):', len(zlib.compress(genome.encode("utf-8"))))
print('CpG Islands:', genome.find('CGGCGG'))

res = translate(genome, codon_table)
# print(res.split('*'))

# Compute protparams
index = 0
for pid, peptide in enumerate(res.split('*')):
    if pid==index:

        n_terminus = lookup_amino(peptide[0])
        c_terminus = lookup_amino(peptide[-1])

        length = len(peptide)

        if length >= 2 and length <= 20:
            type = "Oligopeptide"
        elif length > 20:
            type = "Polypeptide"

        Mw = lookup_weight(peptide)
        net = charge_at_pH(7.0, peptide)
        pI = isoelectric_point(peptide)
        hl = lookup_halflife(peptide[0])
        formula, nb_atoms = atomic_composition(peptide)
        aa_content = amino_count(peptide)
        pos, neg = charged_residues(peptide)
        ec = extinction_coefficient(peptide)
        II = instability_index(peptide)
        ai = aliphatic_index(peptide)
        hp = hydropathy_index(peptide)

        print("Chain Search:", res.find(peptide))
        print(peptide)

        print("N-Terminus:", n_terminus, "| C-Terminus:", c_terminus)
        print("Sequence ID:", pid,
              "| Length:", length,
              "| Type:", type,
              "| Molecular Weight (Da):", round(Mw, 2),
              "| Net Charge (pH = 7.0):", round(net, 2),
              "| Theoretical pI:", round(pI, 2),
              "| Half-life (N-end):", hl)

        print("Atomic Formula:", formula, "| Number of Atoms:", nb_atoms)

        print(aa_content)
        for a, c in aa_content.items():
            print(a, round(c * (100.0/length), 1), "%")

        print("+ charged residues (Arg | Lys | His):", charged_residues(peptide)[0])
        print("- charged residues (Asp | Glu | Cys | Tyr):", charged_residues(peptide)[1])

        if ec == 0:
            print("As there are no Trp, Tyr or Cys in the region considered, this protein should not be visible by UV spectrophotometry.")
        else:
            if peptide.find('W') == -1:
                print("This protein does not contain any Trp residues. Experience shows that this could result in more than 10% error in the computed extinction coefficient.")

            print("Extinction coefficients are in units of M-1 cm-1, at 280nm measured in water.")
            print("Ext. coefficient:", ec)
            print("Abs 0.1% (=1 g/l):", round(ec/Mw, 3))

        print("Instability Index (II):", round(II, 3))
        if II <= 40:
            print("This classifies the protein as stable.")
        else:
            print("This classifies the protein as unstable.")

        print("Aliphatic Index:", round(ai, 3))
        print("Hydropathicity Index (GRAND Average):", round(hp, 3))
