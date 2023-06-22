#!/usr/bin/python

from Bio import SeqIO
import matplotlib.pyplot as plt
def fullsequences(referencefile, alignedfile):

    list1=[]
    ratios =[]
    a =SeqIO.parse(alignedfile, 'fasta')
    ref=SeqIO.parse(referencefile,"genbank")
    refseq = str(ref.__next__().seq)
    records =[]
    b =0

    for i in a:

        b +=1
        string = (str(i.seq))
        newstring = (string.replace('N',''))
        ratio =len(newstring)/len(refseq)
        ratios.append(ratio)

        if len(newstring)/len(refseq)>= 0.95: list1.append(i)
    for i in range(1,b+1): records.append(i)

    ratios.sort(reverse=True)
    
    plt.barh(records,ratios)
    plt.title('Sequence Coverage of Genome')
    plt.ylabel('sequences')
    plt.xlabel('genome covered')
    plt.savefig('Sequence Coverage of Genome.png')
    plt.show()
    SeqIO.write(list1, "longsequencesaligned.fasta","fasta")
print(fullsequences('areference.gb','alignedandsorted.fasta'))