#!/usr/bin/env python3

import yaml, requests

seqs = yaml.safe_load(requests.get("https://www.ncbi.nlm.nih.gov/core/assets/genbank/files/ncov-sequences.yaml").text)
