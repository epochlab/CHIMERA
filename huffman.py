#!/usr/bin/env python3

from libtools import *

UID = 'NC_001542.1'
label, genome = load('genome/' + UID + '.fasta')

def frequency_sort(genome):
    freq = {}
    for n in genome:
        if n in freq:
            freq[n] += 1
        else:
            freq[n] = 1

    chars = freq.keys()
    tuples = []
    for n in chars:
        tuples.append((freq[n], n))
    tuples.sort()
    return tuples

def buildTree(tuples):
    while len(tuples)>1:
        x0 = tuple(tuples[0:2])
        x1 = tuples[2:]
        x = x0[0][0] + x0[1][0]
        tuples = x1 + [(x, x0)]
        tuples.sort()
    return tuples[0]

def trimTree (tree) :
    p = tree[1]
    if type(p) == type(""):
        return p
    else:
        return(trimTree(p[0]), trimTree(p[1]))

queue = frequency_sort(genome)
tree = buildTree(queue)
trim = trimTree(tree)

print(tree)
print(trim)
