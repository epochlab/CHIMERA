#!/usr/bin/env python3

import requests, argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('-uid', type=str, default='NC_045512.2')
args = parser.parse_args(sys.argv[1:])

UID = args.uid

content = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={fasta}&rettype=fasta".format(fasta=UID))
content.raise_for_status()

filename = "genome/" + UID + ".fasta"

with open(filename, 'w') as f:
    f.write(content.text)

print("Download complete", UID)
