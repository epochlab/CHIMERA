#!/usr/bin/env python3

import string

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
