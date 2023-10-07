#!/usr/bin/env python3

import numpy as np

from libtools import *
from molecule import mRNA_codon

UID_1 = 'NC_001542.1'
UID_2 = 'MN996532.2'

_, genome1 = load('genome/' + UID_1 + '.fasta')
_, genome2 = load('genome/' + UID_2 + '.fasta')

n = len(genome1)
m = len(genome2)

score_matrix = np.zeros((n+1, m+1))
traceback_matrix = np.zeros((n+1, m+1))

# Initialize first row and column of scoring matrix
score_matrix[0, 1:] = np.arange(-1, -m-1, -1)
score_matrix[1:, 0] = np.arange(-1, -n-1, -1)

# Fill in scoring matrix and traceback matrix using NumPy array operations
for i in range(1, n+1):
    for j in range(1, m+1):
        match = score_matrix[i-1, j-1] + (genome1[i-1] == genome2[j-1])
        delete = score_matrix[i-1, j] - 1
        insert = score_matrix[i, j-1] - 1
        scores = np.array([match, delete, insert])
        score_matrix[i, j] = np.max(scores)
        traceback_matrix[i, j] = np.argmax(scores)

# Trace back through the traceback matrix to compute alignment
i = n
j = m
alignment1 = ""
alignment2 = ""
while i > 0 or j > 0:
    if traceback_matrix[i, j] == 0:
        alignment1 = genome1[i-1] + alignment1
        alignment2 = genome2[j-1] + alignment2
        i -= 1
        j -= 1
    elif traceback_matrix[i, j] == 1:
        alignment1 = genome1[i-1] + alignment1
        alignment2 = "-" + alignment2
        i -= 1
    else:
        alignment1 = "-" + alignment1
        alignment2 = genome2[j-1] + alignment2
        j -= 1

# Compute similarity score
matches = sum(1 for a, b in zip(alignment1, alignment2) if a == b)
similarity = matches / len(genome1) * 100

# Print out similarity score and aligned sequences
print("Similarity score: {:.2f}%".format(similarity))
print(alignment1)
print(alignment2)