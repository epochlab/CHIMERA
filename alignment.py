#!/usr/bin/env python3

import numpy as np
from tqdm import tqdm
from chimera.fasta import FastaIO

f1 = FastaIO('NC_005831.2')
f2 = FastaIO('NC_006213.1')

n = len(f1.genome)
m = len(f2.genome)

score_matrix = np.zeros((n+1, m+1))
traceback_matrix = np.zeros((n+1, m+1))

# Gap penalties
score_matrix[0, 1:] = np.arange(-1, -m-1, -1)
score_matrix[1:, 0] = np.arange(-1, -n-1, -1)

# Needleman-Wunsch
for i in tqdm(range(1, n+1)):
    for j in range(1, m+1):
        match = score_matrix[i-1, j-1] + (f1.genome[i-1] == f2.genome[j-1])
        delete = score_matrix[i-1, j] - 1
        insert = score_matrix[i, j-1] - 1
        scores = np.array([match, delete, insert])
        score_matrix[i, j] = np.max(scores)
        traceback_matrix[i, j] = np.argmax(scores)

# Traceback to reconstruct aligned sequences (optimal path)
i = n
j = m
alignment1 = ""
alignment2 = ""
while i > 0 or j > 0:
    if traceback_matrix[i, j] == 0:
        alignment1 = f1.genome[i-1] + alignment1
        alignment2 = f2.genome[j-1] + alignment2
        i -= 1
        j -= 1
    elif traceback_matrix[i, j] == 1:
        alignment1 = f1.genome[i-1] + alignment1
        alignment2 = "-" + alignment2
        i -= 1
    else:
        alignment1 = "-" + alignment1
        alignment2 = f2.genome[j-1] + alignment2
        j -= 1

# Compute similarity score
matches = sum(1 for a, b in zip(alignment1, alignment2) if a == b)
similarity = matches / len(f1.genome) * 100

print(f"Similarity score: {similarity:.2f}%")
print("\n" + alignment1)
print("\n" + alignment2)