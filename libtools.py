#!/usr/bin/env python3

import zlib

def compress(genome):
    return zlib.compress(genome.encode("utf-8"))

def translate(seq, dict):
    polypeptide = ""
    for i in range(0, len(seq)-3, 3):
        codon = seq[i:i + 3].replace('T', 'U')
        amino = [k for k, v in dict.items() if codon in v]
        char = str(amino).split('/')[1].replace("']", "").strip()
        polypeptide += str(char)
    return polypeptide
