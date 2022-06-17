#!/usr/bin/env python3

import string

# def binary_encoding(seq):
#     encoded = ''
#     for i in seq:
#         if i == "A":
#             encoded += "00"
#         if i == "G":
#             encoded += "01"
#         if i == "C":
#             encoded += "10"
#         if i == "T":
#             encoded += "11"
#     return encoded

chars = string.ascii_uppercase + string.digits

binary_str = ''.join(format(x, '08b') for x in bytearray(chars, 'utf-8'))
binary_list = [binary_str[i: i+2] for i in range(0, len(binary_str), 2)]

DNA_encoding = {
    "00": "A",
    "01": "G",
    "10": "C",
    "11": "T"
}

nucleo = []
for num in binary_list:
    for key in list(DNA_encoding.keys()):
        if num == key:
            nucleo.append(DNA_encoding.get(key))

DNA = "".join(nucleo)

print(chars)
print(DNA)
