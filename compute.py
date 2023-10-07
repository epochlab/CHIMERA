#!/usr/bin/env python3

from chimera.fasta import FastaIO
from chimera.measure import Molecule, Signal
from chimera.codon import halflife

def main(UID):
    FASTA = FastaIO(UID)

    print("\n" + FASTA.label.upper())
    print("\n" + f"Nucleobases: {len(FASTA.genome)}")
    print(f"GC-Content: {Molecule().gc_content(FASTA.genome):.3f} %")
    print(f"Compression (zlib): {Signal().compress(FASTA.genome)}")
    print(f"Average Hash: {Signal().average_hash(Signal().seq_to_pixels(FASTA.genome))}")

    print("\n" + ">> GENOME PROFILE")
    print(FASTA.genome)

    print("\n" + f">> RESIDUE CHAIN {len(FASTA.res)}")
    f_res = list(filter(None, FASTA.res.split('*')))
    print(f_res)

    index = 0
    for pid, peptide in enumerate(f_res):
        if pid==index:
            n_terminus = Molecule().lookup_amino(peptide[1])
            c_terminus = Molecule().lookup_amino(peptide[-1])
            Mw = Molecule().molecular_weight(peptide)
            net = Molecule().charge_at_pH(7.0, peptide)
            pI = Molecule().isoelectric_point(peptide)
            hl = Molecule().lookup_value(peptide[0], halflife()) # Needs to address whole FASTA
            formula, nb_atoms = Molecule().atomic_composition(peptide)
            aa_content = Molecule().amino_count(peptide)
            ec = Molecule().extinction_coefficient(peptide)
            II = Molecule().instability_index(peptide)
            ai = Molecule().aliphatic_index(peptide)
            hp = Molecule().hydropathy_index(peptide)

            print("\n" + ">> PEPTIDE ANALSIS")
            print("Chain Search:", FASTA.res.find(peptide))
            print(peptide)
            print(f"Sequence ID: {pid}",
                f"| Length: {len(peptide)}",
                f"| Type: {Molecule().peptype(peptide)}",
                f"| Molecular Weight (Da): {Mw:.2f}",
                f"| Net Charge (pH = 7.0): {net:.2f}",
                f"| Theoretical pI: {pI:.2f}",
                f"| Half-life (N-end): {hl}")

            print(f"N-Terminus: {n_terminus} | C-Terminus: {c_terminus}")
            print(f"Atomic Formula: {formula}| Number of Atoms: {nb_atoms}")

            for a, c in aa_content.items():
                print(f"{a} {c * (100.0/len(peptide)):.1f} % {c}")

            print(f"+ charged residues (Arg | Lys | His): {Molecule().charged_residues(peptide)[0]}")
            print(f"- charged residues (Asp | Glu | Cys | Tyr): {Molecule().charged_residues(peptide)[1]}")

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

if __name__ == "__main__":
    main('NC_001542.1') # RABIES VIRUS, COMPLETE GENOME