#!/usr/bin/env python3

import requests

UID = "NC_045512.2"
content = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={fasta}&rettype=fasta".format(fasta=UID))
content.raise_for_status()

filename = "data/" + UID + ".txt"

with open(filename, 'w') as f:
    f.write(content.text)
    print("Download complete", UID)
