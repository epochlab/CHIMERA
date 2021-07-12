#!/usr/bin/env python3
import zlib

nucleotides = ['A', 'C', 'G', 'T']

file = open('genome.txt')

data = str(file.readlines()[692:1191])
for s in "[]\ n,'0123456789":
    data = data.replace(s, '')
data = data.upper()

print(data)

length = len(data)
dst = zlib.compress(data.encode("utf-8"))

print('Nucleotides:', length)
print('Compression:', len(dst))
