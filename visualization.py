#!/usr/bin/env python3

from libtools import *
from dict import mRNA_codon

codon_table = mRNA_codon()

fasta = 'NC_001563.2'
label, genome = load("genome/" + fasta + ".txt")
binary = binary_encoding(genome)

print(binary)
