#!/usr/bin/python

from Bio import SeqIO
import numpy as np

with open('alignedandsorted.fasta') as f:
    file =SeqIO.parse(f, 'fasta')
    listofL=[]
    startL = 8568
    endL = 15065
    for record in file:
        listofL.append(record[startL:endL])

    SeqIO.write(listofL, 'onlySN1-M%GENE.fasta', 'fasta')


