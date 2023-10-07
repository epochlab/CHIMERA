#!/usr/bin/env python3

from libtools import FastaIO, Molecule, Signal
from codon import RNA, halflife

FASTA = FastaIO()
MOLECULE = Molecule()
SIGNAL = Signal()

def main(UID):
    label, genome = FASTA.load(f'genome/{UID}.fasta')

    res = FASTA.translate(genome, RNA())
    f_res = list(filter(None,res.split('*')))

    print("\n" + label.upper())
    print("\n" + f"Nucleobases: {len(genome)}")
    print(f"GC-Content: {MOLECULE.gc_content(genome):.3f} %")
    print(f"Compression (zlib): {SIGNAL.compress(genome)}")
    print(f"Average Hash: {SIGNAL.average_hash(SIGNAL.seq_to_pixels(genome))}")

    print("\n" + ">> GENOME PROFILE")
    print(genome)

    print("\n" + f">> RESIDUE CHAIN {len(res)}")
    print(f_res)

    index = 0
    for pid, peptide in enumerate(f_res):
        if pid==index:
            n_terminus = MOLECULE.lookup_amino(peptide[1])
            c_terminus = MOLECULE.lookup_amino(peptide[-1])
            Mw = MOLECULE.molecular_weight(peptide)
            net = MOLECULE.charge_at_pH(7.0, peptide)
            pI = MOLECULE.isoelectric_point(peptide)
            hl = MOLECULE.lookup_value(peptide[0], halflife()) # Needs to address whole FASTA
            formula, nb_atoms = MOLECULE.atomic_composition(peptide)
            aa_content = MOLECULE.amino_count(peptide)
            ec = MOLECULE.extinction_coefficient(peptide)
            II = MOLECULE.instability_index(peptide)
            ai = MOLECULE.aliphatic_index(peptide)
            hp = MOLECULE.hydropathy_index(peptide)

            print("\n" + ">> PEPTIDE ANALSIS")
            print("Chain Search:", res.find(peptide))
            print(peptide)
            print(f"Sequence ID: {pid}",
                f"| Length: {len(peptide)}",
                f"| Type: {MOLECULE.peptype(peptide)}",
                f"| Molecular Weight (Da): {Mw:.2f}",
                f"| Net Charge (pH = 7.0): {net:.2f}",
                f"| Theoretical pI: {pI:.2f}",
                f"| Half-life (N-end): {hl}")

            print(f"N-Terminus: {n_terminus} | C-Terminus: {c_terminus}")
            print(f"Atomic Formula: {formula}| Number of Atoms: {nb_atoms}")

            for a, c in aa_content.items():
                print(f"{a} {c * (100.0/len(peptide)):.1f} % {c}")

            print(f"+ charged residues (Arg | Lys | His): {MOLECULE.charged_residues(peptide)[0]}")
            print(f"- charged residues (Asp | Glu | Cys | Tyr): {MOLECULE.charged_residues(peptide)[1]}")

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

UID = 'NC_001542.1' # RABIES VIRUS, COMPLETE GENOME
if __name__ == "__main__":
    main(UID)