#!/usr/bin/python

from Bio import SeqIO
import pandas as pd
seq = SeqIO.parse("data/concatenated.fasta","fasta")
cleaned_seq =[]
data ={}
for record in seq:
    name, acc, date = record.id.split('|')
    record.id = name
    record.description = ""
    cleaned_seq.append(record)
    virus, subtype, country, id, year =name.split('/')
    data[name]={"strain": name, "virus": virus,"subtype": subtype, "country":country, "acc":acc,"date":date,}
with open("data/sequences.fasta","w") as f:
        SeqIO.write(cleaned_seq, f, "fasta")
df = pd.DataFrame.from_dict(data, orient='index')
df.to_csv("data/metadata.tsv",sep="\t",index=False)
