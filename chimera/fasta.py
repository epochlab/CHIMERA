#!/usr/bin/env python3

from chimera.codon import RNA

class FastaIO():
    def __init__(self, UID):
        self.UID = UID
        self.path = f"genome/{self.UID}.fasta"
        self.codon = RNA()
        
        self.label, self.genome = self.load(self.path)
        self.res = self.translate(self.genome)

    def load(self, fasta):
        with open(fasta) as f:
            header = f.readline().rstrip()
            seq = ''.join(line.strip() for line in f)
        return header, seq

    def transcribe(self, seq):
        return seq.replace('T', 'U') # DNA > RNA transcription - Thymine (T) is replaced with Uracil (U).

    def translate(self, seq):
        i, count = 0, 1
        res = ''
        while i < len(seq):
            codon = self.transcribe(seq[i:i+3])
            amino = [k for k, v in self.codon.items() if codon in v]

            if codon=='AUG':                                                        # START open reading frame
                count = 3
            if count==3 and len(codon)==3:
                res += str(amino).split('/')[1].replace("']", "").strip()
            if codon=='UAG' or codon=='UAA' or codon=='UGA':                        # STOP open reading frame
                i += 2
                count = 1

            i += count
        return res