#!/usr/bin/python


import json
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
def KnownGenomePlot (filename):
    with open('aa_muts.json') as f:
        aamuts = json.load(f)
        dictofgenes=dict()
        keys = ('NS1','NS2','N','P','M','SH','G','F','M2','L')

        for i in keys:
            for key, node in aamuts['annotations'].items():
                if key == i:
                    g=[]
                    g = list(range(node['start'],node['end']+1))
            dictofgenes[i]=g
    n = 0
    x = []
    mylist1 = []
    mylist2 =[]
    records = SeqIO.parse(filename, "fasta")
    align_array = np.array([record.seq for record in records])
    for line in align_array:
        n+=1
    a = (align_array.T)
    for element in a: 
        unique, counts = np.unique(element, return_counts=True)
        b =dict(zip(unique, counts))
        c = b.get('N')
        mylist1.append(c)
    for i in mylist1:
        d = n-i
        mylist2.append(d)

    fig = plt.figure()
    for i in range(len(mylist2)):
        x.append(i)
    height = mylist2
    width = 1
    colourslist=[]
    colourlist =['tomato', 'yellow', 'black', 'deepskyblue', 'bisque', 'violet', 'forestgreen', 'indigo', 'cyan','darkred']
    for i in x: colourslist.append('gray')
    for (gene,number), colour in zip(dictofgenes.items(), colourlist):
        for i in x:
            if i in number:
                colourslist[i-1]=colour
    patches =[]
    for colour, gene in zip(colourlist, keys):
        patch = mpatches.Patch(color=colour, label=gene)
        patches.append(patch)
    plt.legend(handles=patches)
    plt.bar(x, height, width, color=colourslist)
    plt.title("Known Sequences")
    plt.xlabel('Genome')
    plt.ylabel('Sequence Number')
    plt.savefig('KnownSequences.pdf')

print(KnownGenomePlot("alignedandsorted.fasta"))
