#!/usr/bin/env python3

from PIL import Image
import numpy as np
import math

from libtools import *

def binary_array_to_hex(arr):
	bit_string = ''.join(str(b) for b in 1 * arr.flatten())
	width = int(np.ceil(len(bit_string)/4))
	return '{:0>{width}x}'.format(int(bit_string, 2), width=width)

def average_hash(image, hash_size=8, mean=np.mean):
    image = image.convert('L').resize((hash_size, hash_size), Image.ANTIALIAS)
    pixels = np.asarray(image)
    diff = pixels > mean(pixels)
    return binary_array_to_hex(diff)

UID = 'NC_001474.2'
label, genome = load('genome/' + UID + '.fasta')

pixels = []

for n in genome:
    if n == 'A':
        pixels.append((10,10,10))
    if n == 'C':
        pixels.append((255,255,255))
    if n == 'G':
        pixels.append((65,150,65))
    if n == 'T':
        pixels.append((0,50,140))

array = np.array(pixels, dtype=np.uint8)

width = round(math.sqrt(len(genome)))
height = len(genome) // width

img = Image.frombytes("RGB", (width, height), bytes(array))

hash = average_hash(img)

img.save(UID + "_" + hash + '.png')

print("Hashing complete:", hash)
