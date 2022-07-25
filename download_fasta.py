#!/usr/bin/env python3

import requests, argparse, sys, os, csv
import pandas as pd
from libtools import *

parser = argparse.ArgumentParser()
parser.add_argument('-uid', type=str, default='NC_001542.1')
args = parser.parse_args(sys.argv[1:])

UID = args.uid

database = 'genome_db.csv'
valid = os.path.exists(database)

fieldnames = ['uid', 'name', 'length', 'zlib', 'hash']

reader = pd.read_csv(database)
# print(reader.to_string())

stored = False
for row in reader["uid"]:
    if row == UID:
        stored = True
        print(row, "found in library.")
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

    new_row = pd.DataFrame(
        {'uid': UID,
        'name': (" ").join(label.split(" ")[1:]),
        'length': length,
        'zlib': size,
        'hash': hash}, index=[0]
    )

    new_frame = pd.concat([new_row, reader.loc[:]])
    new_frame = new_frame.sort_values(by='uid').reset_index(drop=True)
    new_frame.to_csv(database, index=False)

    print(label)
    print(length, size, hash)

    # print(new_frame.to_string())