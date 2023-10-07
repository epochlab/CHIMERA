#!/usr/bin/env python3

import requests, argparse, sys, os
import pandas as pd
from libtools import FastaIO, Signal

parser = argparse.ArgumentParser()
parser.add_argument('-uid', type=str, default='NC_001542.1')
args = parser.parse_args(sys.argv[1:])

UID = args.uid

db = 'genome_db.csv'

# Check database exists
if os.path.exists(db):
    reader = pd.read_csv(db)
else:
    print("No directory found")

# Check for duplicates
for row in reader["uid"]:
    if row == UID:
        stored = True
        print(row, "found in library")
    else:
        stored = False

# Download and store FASTA
if stored == False:
    content = requests.get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={UID}&rettype=fasta")
    content.raise_for_status()

    filename = f"genome/{UID}.fasta"
    with open(filename, 'w') as f:
        f.write(content.text)

    label, genome = FastaIO().load(filename)

    pixels = Signal().seq_to_pixels(genome)
    hash =  Signal().average_hash(pixels)

    new_row = pd.DataFrame(
        {'uid': UID,
        'name': (" ").join(label.split(" ")[1:]),
        'length': len(genome),
        'zlib': Signal().compress(genome),
        'hash': hash}, index=[0]
        )

    new_frame = pd.concat([new_row, reader.loc[:]])
    new_frame = new_frame.sort_values(by='uid').reset_index(drop=True)
    
    new_frame.to_csv(db, index=False)
    print(f"{label} added to library")