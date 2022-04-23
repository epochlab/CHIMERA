#!/usr/bin/env python3

import os, csv
from libtools import *

UID = 'NC_045512.2'
label, genome = load('genome/' + UID + '.fasta')

pixels = seq_to_pixels(genome)

length = len(genome)
size = compress(genome)
hash =  average_hash(pixels)

fieldnames = ['uid', 'name', 'length', 'size', 'hash']

rows = [
    {'uid': UID,
    'name': (" ").join(label.split(" ")[1:]),
    'length': length,
    'size': size,
    'hash': hash},
]

database = 'hash_db.csv'
valid = os.path.exists(database)

with open(database, 'a', encoding='UTF8') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    if  valid == False:
        writer.writeheader()
    writer.writerows(rows)

print(label)
print(length, size, hash)

# pixels.save(UID + "_" + hash + '.png')
