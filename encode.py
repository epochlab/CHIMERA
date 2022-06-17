#!/usr/bin/env python3

import string

chars = string.ascii_uppercase + string.digits

binary_str = ''.join(format(x, '08b') for x in bytearray(chars, 'utf-8'))
binary_list = [binary_str[i: i+2] for i in range(0, len(binary_str), 2)]

dict = {
    "00": "A",
    "01": "G",
    "10": "C",
    "11": "T"
    }

string = []
for n in binary_list:
    for key in list(dict.keys()):
        if n == key:
            string.append(dict.get(key))

DNA = "".join(string)

print(chars)
print(binary_str)
print(DNA)
