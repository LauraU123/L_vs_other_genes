#!/usr/bin/python

import numpy as np
from Bio import SeqIO

def sorted_file(referencefile, alignedfile):
    refseq = SeqIO.parse(referencefile, "genbank")
    seq2 = np.array(refseq.__next__())
    list1 =[]
    for record in SeqIO.parse(alignedfile,"fasta"):

        seq1 = np.array(record.seq)
        good_indices = (refseq!='N')&(seq1!='N')
        a = np.mean(seq1[good_indices]==seq2[good_indices])

        if a >0.9: list1.append(record)
        else: continue

    SeqIO.write(list1, "alignedandsorted.fasta","fasta")
    print(len(list1))

print(sorted_file('areference.gb','aligned.fasta'))
