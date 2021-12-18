#!/usr/bin/env python3

      # A: Ala, name: Alanine
      # I: Ile, name: Isoleucine
      # R: Arg, name: Arginine
      # L: Leu, name: Leucine
      # N: Asn, name: Asparagine
      # K: Lys, name: Lysine
      # D: Asp, name: Aspartic Acid
      # M: Met, name: Methionine
      # F: Phe, name: Phenylalanine
      # C: Cys, name: Cysteine
      # P: Pro, name: Proline
      # Q: Gln, name: Glutamine
      # S: Ser, name: Serine
      # E: Glu, name: Glutamic Acid
      # T: Thr, name: Threonine
      # W: Trp, name: Tryptophan
      # G: Gly, name: Glycine
      # Y: Tyr, name: Tyrosine
      # H: His, name: Histidine
      # V: Val, name: Valine
      #
      # U: Sec, name: Selenocysteine # uga (stop)
      # O: Pyl  # uag (stop)
      #
      # B: Asx  # D/N
      # Z: Glx  # E/Q
      # J: Xle  # I/L # => Illigal in protparams
      # X: Xaa  # unknown

nucleotides = ['A', 'C', 'G', 'T']

def mRNA_codon():
    dict = {
        'Ala / A': ['GCU', 'GCC', 'GCA', 'GCG'],
        'Ile / I': ['AUU', 'AUC', 'AUA'],
        'Arg / R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'Leu / L': ['CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG'],
        'Asn / N': ['AAU', 'AAC'],
        'Lys / K': ['AAA', 'AAG'],
        'Asp / D': ['GAU', 'GAC'],
        'Met / M': ['AUG'],
        'Phe / F': ['UUU', 'UUC'],
        'Cys / C': ['UGU', 'UGC'],
        'Pro / P': ['CCU', 'CCC', 'CCA', 'CCG'],
        'Gln / Q': ['CAA', 'CAG'],
        'Ser / S': ['UCU', 'UCC', 'UCA', 'UCG',  'AGU', 'AGC'],
        'Glu / E': ['GAA', 'GAG'],
        'Thr / T': ['ACU', 'ACC', 'ACA', 'ACG'],
        'Trp / W': ['UGG'],
        'Gly / G': ['GGU', 'GGC', 'GGA', 'GGG'],
        'Tyr / Y': ['UAU', 'UAC'],
        'His / H': ['CAU', 'CAC'],
        'Val / V': ['GUU', 'GUC', 'GUA', 'GUG'],
        'STOP / *': ['UAA', 'UGA', 'UAG'],
        }
    return dict

def molecular_weight():
    dict = {
        'Ala / A': [89.09],
        'Ile / I': [131.18],
        'Arg / R': [174.20],
        'Leu / L': [131.18],
        'Asn / N': [132.12],
        'Lys / K': [146.19],
        'Asp / D': [133.10],
        'Met / M': [149.21],
        'Phe / F': [165.19],
        'Cys / C': [121.16],
        'Pro / P': [115.13],
        'Gln / Q': [146.15],
        'Ser / S': [105.09],
        'Glu / E': [147.13],
        'Thr / T': [119.12],
        'Trp / W': [204.23],
        'Gly / G': [75.07],
        'Tyr / Y': [181.19],
        'His / H': [155.16],
        'Val / V': [117.15],
        'STOP / *': [],
        }
    return dict

# HostID: Mammalian | Yeast | E.Coli
# Unit: Minutes
def halflife():
    dict = {
        'Ala / A': [264, 1200, 600],
        'Ile / I': [1200, 30, 600],
        'Arg / R': [60, 2, 2],
        'Leu / L': [330, 3, 2],
        'Asn / N': [84, 3, 600],
        'Lys / K': [78, 3, 2],
        'Asp / D': [66, 3, 600],
        'Met / M': [1800, 1200, 600],
        'Phe / F': [66, 3, 2],
        'Cys / C': [72, 1200, 600],
        'Pro / P': [1200, 1200, 600],
        'Gln / Q': [48, 10, 600],
        'Ser / S': [114, 1200, 600],
        'Glu / E': [60, 30, 600],
        'Thr / T': [432, 1200, 600],
        'Trp / W': [168, 3, 2],
        'Gly / G': [1800, 1200, 600],
        'Tyr / Y': [168, 10, 2],
        'His / H': [210, 10, 600],
        'Val / V': [6000, 1200, 600],
        'STOP / *': [],
        }
    return dict

def hydropathy():
    dict = {
        'Ala / A': [1.8],
        'Ile / I': [4.5],
        'Arg / R': [-4.5],
        'Leu / L': [3.8],
        'Asn / N': [-3.5],
        'Lys / K': [-3.9],
        'Asp / D': [-3.5],
        'Met / M': [1.9],
        'Phe / F': [2.8],
        'Cys / C': [2.5],
        'Pro / P': [-1.6],
        'Gln / Q': [-3.5],
        'Ser / S': [-0.8],
        'Glu / E': [-3.5],
        'Thr / T': [-0.7],
        'Trp / W': [-0.9],
        'Gly / G': [-0.4],
        'Tyr / Y': [-1.3],
        'His / H': [-3.2],
        'Val / V': [4.2],
        'STOP / *': [],
        }
    return dict

def intrinsic_disorder():
    dict = {
        'Ala / A': [0.06],
        'Ile / I': [-0.486],
        'Arg / R': [0.180],
        'Leu / L': [-0.326],
        'Asn / N': [0.007],
        'Lys / K': [0.586],
        'Asp / D': [0.192],
        'Met / M': [-0.397],
        'Phe / F': [-0.697],
        'Cys / C': [0.02],
        'Pro / P': [0.987],
        'Gln / Q': [0.318],
        'Ser / S': [0.341],
        'Glu / E': [0.736],
        'Thr / T': [0.059],
        'Trp / W': [-0.884],
        'Gly / G': [0.166],
        'Tyr / Y': [-0.510],
        'His / H': [0.303],
        'Val / V': [-0.121],
        'STOP / *': [0.02],
        }
    return dict

def average_mass():
    dict = {
        'Ala / A': [71.0788],
        'Ile / I': [113.1594],
        'Arg / R': [156.1875],
        'Leu / L': [113.1594],
        'Asn / N': [114.1038],
        'Lys / K': [128.1741],
        'Asp / D': [115.0886],
        'Met / M': [131.1926],
        'Phe / F': [147.1766],
        'Cys / C': [103.1388],
        'Pro / P': [97.1167],
        'Gln / Q': [128.1307],
        'Ser / S': [87.0782],
        'Glu / E': [129.1155],
        'Thr / T': [101.1051],
        'Trp / W': [186.2132],
        'Gly / G': [57.0519],
        'Tyr / Y': [163.1760],
        'His / H': [137.1411],
        'Val / V': [99.1326],
        'STOP / *': [150.0388],
        }
    return dict

def atom():    # C = Carbon | H = Hydrogen | N = Nitrogen | O = Oxygen | S = Sulfur
    dict = {
        'Ala / A': (3, 7, 1, 2, 0), # C3H7NO2
        'Ile / I': [6, 13, 1, 2, 0], # C6H13NO2
        'Arg / R': [6, 14, 4, 2, 0], # C6H14N4O2
        'Leu / L': [6, 13, 1, 2, 0], # C6H13NO2
        'Asn / N': [4, 8, 2, 3, 0], # C4H8N2O3
        'Lys / K': [6, 14, 2, 2, 0], # C6H14N2O2
        'Asp / D': [4, 7, 1, 4, 0], # C4H7NO4
        'Met / M': [5, 11, 1, 2, 1], # C5H11NO2S
        'Phe / F': [9, 11, 1, 2, 0], # C9H11NO2
        'Cys / C': [3, 7, 1, 2, 1], # C3H7NO2S
        'Pro / P': [5, 9, 1, 2, 0], # C5H9NO2
        'Gln / Q': [5, 10, 2, 3, 0], # C5H10N2O3
        'Ser / S': [3, 7, 1, 3, 0], # C3H7NO3
        'Glu / E': [5, 9, 1, 4, 0], # C5H9NO4
        'Thr / T': [4, 9, 1, 3, 0], # C4H9NO3
        'Trp / W': [11, 12, 2, 2, 0], # C11H12N2O2
        'Gly / G': [2, 5, 1, 2, 0], # C2H5NO2
        'Tyr / Y': [9, 11, 1, 3, 0], # C9H11NO3
        'His / H': [6, 9, 3, 2, 0], # C6H9N3O2
        'Val / V': [5, 11, 1, 2, 0], # C5H11NO2
        'STOP / *': [],
        }
    return dict
