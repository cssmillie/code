#!/usr/bin/env python

import sys, util

codon_table = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
               "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
               "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
               "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
               "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
               "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
               "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
               "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
               "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
               "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
               "TAA" : "X", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
               "TAG" : "X", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
               "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
               "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
               "TGA" : "X", "CGA" : "R", "AGA" : "R", "GGA" : "G",
               "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
               }


def get_codons(nt):
    if len(nt) % 3 != 0:
        sys.stderr.write('error: len(nt) not divisible by 3\n')
        return False
    for i in range(int(len(nt)/3)):
        beg = i*3
        end = beg + 3
        yield nt[beg:end]

def translate(nt):
    nt = nt.upper()
    aa = ''
    for codon in get_codons(nt):
        if codon is False:
            return False
        if codon in codon_table:
            aa += codon_table[codon]
        else:
            aa += 'X'
    return aa

if __name__ == '__main__':
    import sys
    x = {}
    if len(sys.argv) > 1:
        fh = sys.argv[1]
    else:
        fh = sys.stdin
    for record in util.iter_fst(fh):
        if record[0] in x:
            continue
        x[record[0]] = 1
        record[1] = translate(record[1])
        if record[1]:
            print('\n'.join(record))
