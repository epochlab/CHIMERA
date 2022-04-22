#!/usr/bin/env python3

from PIL import Image
import numpy as np
import math

from libtools import *

UID = 'NC_001474.2'
label, genome = load('genome/' + UID + '.fasta')

pixels = []

for n in genome:
    if n == 'A':
        pixels.append((255,0,0))
    if n == 'C':
        pixels.append((0,255,0))
    if n == 'G':
        pixels.append((0,0,255))
    if n == 'T':
        pixels.append((255,255,255))

# print(col)
# print(label, genome)

array = np.array(pixels, dtype=np.uint8)

width = round(math.sqrt(len(pixels)))
height = len(pixels) // width

img = Image.frombytes("RGB", (width, height), bytes(array))

img.save('hash_' + UID + '.png')
