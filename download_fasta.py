#!/usr/bin/env python3

import requests, argparse, sys, os, csv, math
from libtools import *

database = 'hash_db.csv'
valid = os.path.exists(database)

parser = argparse.ArgumentParser()
parser.add_argument('-uid', type=str, default='NC_045512.2')
args = parser.parse_args(sys.argv[1:])

UID = args.uid

stored = False
with open(database, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0] == UID:
            stored = True
            print(row[0], "found in genome library.")
            break

if stored == False:
    content = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={fasta}&rettype=fasta".format(fasta=UID))
    content.raise_for_status()

    filename = "genome/" + UID + ".fasta"

    with open(filename, 'w') as f:
        f.write(content.text)

    label, genome = load(filename)
    pixels = seq_to_pixels(genome)

    length = len(genome)
    size = compress(genome)
    hash =  average_hash(pixels)

    fieldnames = ['uid', 'name', 'length', 'zlib', 'hash']

    rows = [
        {'uid': UID,
        'name': (" ").join(label.split(" ")[1:]),
        'length': length,
        'zlib': size,
        'hash': hash},
    ]

    with open(database, 'a', encoding='UTF8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if  valid == False:
            writer.writeheader()
        writer.writerows(rows)

    print(label)
    print(length, size, hash)

    # pixels.save(UID + "_" + hash + '.png')
