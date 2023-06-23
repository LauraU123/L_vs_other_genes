#!/usr/bin/python

from Bio import SeqIO
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="make new reference depending on whether the entire genome or only part is to be used for the tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="GenBank file with reference sequences")
    parser.add_argument("--output", required=True, help="GenBank new reference file")
    args = parser.parse_args()

    with open(args.alignment) as f:
        file =SeqIO.parse(f, 'fasta')
        listofL=[]
        startL = 8568
        endL = 15065
        for record in file:
            listofL.append(record[startL:endL])

        SeqIO.write(listofL, args.output, 'fasta')


