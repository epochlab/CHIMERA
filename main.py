#!/usr/bin/env python3

import re

data = 'genome.txt'
nucleobases = ['a', 'c', 'g', 't']

sequence = []
with open(data) as f:
    line = sorted(f.readlines()[692:1191])
    sequence.append(line)

sequence = str(sequence)[:100]
print(sequence)
