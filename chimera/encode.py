#!/usr/bin/env python3

def encode(binary, encoding):
    string = []
    for n in [binary[i: i+2] for i in range(0, len(binary), 2)]:
        for key in list(encoding.keys()):
            if n == key:
                string.append(encoding.get(key))
    return "".join(string)

mapping = {
    "00": "A",
    "01": "G",
    "10": "C",
    "11": "T"
    }

data = "ABC123<>?!@£$%^&*()"
binary = ''.join(format(x, '08b') for x in bytearray(data, 'utf-8'))
DNA = encode(binary, mapping)

print(f"Message: {data}")
print(f"Binary: {binary}")
print(f"Residue: {DNA}")