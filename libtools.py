#!/usr/bin/env python3

def load(input):
    nucleotides = ['A', 'C', 'G', 'T']
    file = open(input[0])
    data = str(file.readlines()[input[1][0]:input[1][1]])
    genome = [x for x in data.upper() if (x in nucleotides)]
    genome = ''.join(genome)
    return genome

def translate(seq, dict):
    polypeptide = ""
    for i in range(0, len(seq)-3, 3):
        codon = seq[i:i + 3].replace('T', 'U')                                  # DNA to RNA transcription > Thymine is replace with Uracil
        amino = [k for k, v in dict.items() if codon in v]
        char = str(amino).split('/')[1].replace("']", "").strip()
        polypeptide += str(char)
    return polypeptide
