#!/usr/bin/env python3

from libtools import *
from molecule import mRNA_codon

def main(UID):
    label, genome = load(f'genome/{UID}.fasta')

    res = translate(genome, mRNA_codon())
    f_res = list(filter(None,res.split('*')))

    print("\n" + label.upper())
    print("\n" + f"Nucleobases: {len(genome)}")
    print(f"GC-Content: {gc_content(genome):.3f} %")
    print(f"Compression (zlib): {compress(genome)}")
    print(f"Average Hash: {average_hash(seq_to_pixels(genome))}")

    print("\n" + ">> GENOME PROFILE")
    print(genome)

    print("\n" + f">> RESIDUE CHAIN {len(res)}")
    print(f_res)

    index = 0
    for pid, peptide in enumerate(f_res):
        if pid==index:

            n_terminus = lookup_amino(peptide[1])
            c_terminus = lookup_amino(peptide[-1])
            Mw = molecular_weight(peptide)
            net = charge_at_pH(7.0, peptide)
            pI = isoelectric_point(peptide)
            hl = lookup_value(peptide[0], halflife())
            formula, nb_atoms = atomic_composition(peptide)
            aa_content = amino_count(peptide)
            pos, neg = charged_residues(peptide)
            ec = extinction_coefficient(peptide)
            II = instability_index(peptide)
            ai = aliphatic_index(peptide)
            hp = hydropathy_index(peptide)

            print("\n" + ">> PEPTIDE ANALSIS")
            print("Chain Search:", res.find(peptide))
            print(peptide)
            print(f"Sequence ID: {pid}",
                f"| Length: {len(peptide)}",
                f"| Type: {peptype(peptide)}",
                f"| Molecular Weight (Da): {Mw:.2f}",
                f"| Net Charge (pH = 7.0): {net:.2f}",
                f"| Theoretical pI: {pI:.2f}",
                f"| Half-life (N-end): {hl}")

            print(f"N-Terminus: {n_terminus} | C-Terminus: {c_terminus}")
            print(f"Atomic Formula: {formula}| Number of Atoms: {nb_atoms}")

            for a, c in aa_content.items():
                print(a, round(c * (100.0/len(peptide)), 1), '%', c)

            print(f"+ charged residues (Arg | Lys | His): {charged_residues(peptide)[0]}")
            print(f"- charged residues (Asp | Glu | Cys | Tyr): {charged_residues(peptide)[1]}")

            if ec == 0:
                print("As there are no Trp, Tyr or Cys in the region considered, this protein should not be visible by UV spectrophotometry.")
            else:
                if peptide.find('W') == -1:
                    print("This protein does not contain any Trp residues. Experience shows that this could result in more than 10% error in the computed extinction coefficient.")

                print("Extinction coefficients are in units of M-1 cm-1, at 280nm measured in water.")
                print(f"Ext. coefficient: {ec}")
                print(f"Abs 0.1% (=1 g/l): {(ec/Mw):.3f}")

            print(f"Instability Index (II): {II:.3f}")
            if II <= 40:
                print("This classifies the protein as stable.")
            else:
                print("This classifies the protein as unstable.")

            print(f"Aliphatic Index: {ai:.3f}")
            print(f"Hydropathicity Index (GRAND Average): {hp:.3f}")

UID = 'NC_001542.1'
if __name__ == "__main__":
    main(UID)