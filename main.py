#!/usr/bin/env python3
import zlib

nucleotides = ['A', 'C', 'G', 'T']

file = open('genome.txt')

data = str(file.readlines()[692:1191])
genome = [ele for ele in data.upper() if (ele in nucleotides)]
genome = ''.join(genome)

length = len(genome)
dst = zlib.compress(genome.encode("utf-8"))

print(genome)
print('Base Pairs:', length)
print('Compression:', len(dst))
