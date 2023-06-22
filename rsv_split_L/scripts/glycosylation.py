#!/usr/bin/python

import argparse, json
from random import sample
import numpy as np
from Bio import Phylo, AlignIO, SeqIO
import re

def glycosylation_count(total_aa_seq, glyc_mask=None):
    if glyc_mask is None:
        glyc_mask = np.ones(len(total_aa_seq), dtype=bool)

    total_aa_seq_masked = "".join([aa if mask else 'X'
                                   for (mask, aa) in zip(glyc_mask, total_aa_seq)])

    return len(re.findall('N[^P][ST][^P]', total_aa_seq_masked))


with open ('onlyLL.fasta') as f:
    file = SeqIO.parse(f, 'fasta')
    for i in file:
        print(i.id, glycosylation_count(i.seq))
    
    
    
tree = "results/tree.nwk"
alignment = 'onlyLL.fasta'



T = Phylo.read(tree, 'newick')

glyc_json = {}
aln = {s.name:str(s.seq) for s in AlignIO.read(alignment, 'fasta')}
root_seq = aln[T.root.name]
root_glyc = glycosylation_count(root_seq)
for n in T.find_clades(order='preorder'):
    glyc_json[n.name] = {'glyc':glycosylation_count(aln[n.name]) - root_glyc}

with open('glyc.json', 'wt') as fh:
    json.dump({'nodes':glyc_json}, fh)
