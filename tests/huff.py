#!/usr/bin/env python3

import sys
sys.path.append('../CHIMERA')

from chimera.fasta import FastaIO
from chimera.huffman import *
from collections import Counter

FASTA = FastaIO('NC_001542.1')

print(FASTA.label.upper())
print(f"Nucleobases: {len(FASTA.genome)}")

freq = dict(Counter(FASTA.genome))
count = sorted(freq.items(), key=lambda x: x[1], reverse=True)
print("Count:", count)

node = build_tree(count)
encoding = assign_code(node)

for i in encoding:
    print(f'{i} : {encoding[i]}')

print(encode(FASTA.genome, encoding))