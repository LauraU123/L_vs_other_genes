
import json
import math
import pandas as pd
from collections import Counter

def GenesAndCodons(aamutants, ntmutants):

    with open(aamutants) as f:
        with open(ntmutants) as g:
            list1=[]
            list2=[] 
            list3=[]
            dictofgenes=dict()
            aamuts = json.load(f)
            ntmuts = json.load(g)
            keys = ('F','G','M','M2','NS1','NS2','P','SH','N', 'L')

            for i in keys:
                for key, node in aamuts['annotations'].items():
                    if key == i:
                        g=[]
                        g = list(range(node['start'],node['end']+1))
                dictofgenes[i]=g
            for k, n in aamuts['nodes'].items():
                for key, node in ntmuts['nodes'].items():
                    if k == key:
                        for y in node['muts']:
                            numbers =[]
                            number =int(y[1:-1])
                            numbers.append(number)
                            for x in numbers:
                                for a,g in dictofgenes.items():
                                    if x in g:
                                        list3.append(x)
                                        index = a
                                        codon = (math.floor((x-g[0])/3))+1        
                                        list2.append(index)
                                        list1.append(codon)
        df=pd.DataFrame({'Gene':list2,'Codon':list1, 'Nucleotide from start': list3})
        return(df)

def MutationsineachGene(aamutations, ntmutations):

    genes =['F', 'G', 'M', 'M2', 'NS1', 'NS2', 'P', 'SH', 'N', 'L']
    dictionary = dict()
    dictionary1=dict()
    df=(GenesAndCodons(aamutations, ntmutations))
    for i in genes:

        dictionary[i]= df.loc[df['Gene']==i]

    for i,j in dictionary.items():

        a =(j['Codon'].value_counts())
        b = pd.DataFrame({'Codon':a.index, 'Frequency':a.values})

    for i,j in dictionary.items():
        j =j.reset_index(drop=True)
        dictionary1[i]=j
    return(dictionary1)

def AA_Mutations(aamutations, ntmutations):
    aa_m = dict()
    with open(aamutations) as f:
        with open(ntmutations) as g:

            lis1 =[]
            keys = ('F','G','M','M2','NS1','NS2','P','SH','N', 'L')
            aamuts = json.load(f)

            for a in keys:

                lis1=[]

                for k, n in aamuts['nodes'].items():  
                            for i,j in n['aa_muts'].items():
                                if j!=[] and i ==a: lis1.append(j)
                flatlist =[item for sublist in lis1 for item in sublist]
                flatlist = [int(i[1:-1]) for i in flatlist]
                aa_m[a]=flatlist
    return(aa_m)


mylist=[]
a = AA_Mutations('aa_muts.json', 'nt_muts.json')
nt_muts = MutationsineachGene('aa_muts.json', 'nt_muts.json')['G']
ntmuts_counter =(Counter(list(nt_muts['Codon'])))
for j in a['G']: mylist.append(j)

aa_muts_counter = (Counter(mylist))
synonymous = ntmuts_counter-aa_muts_counter
synonymouslist =[]
for i, j in synonymous.items():
    for n in range(j): synonymouslist.append(i)


with open('aa_muts.json') as f:
    aamuts = json.load(f)

    for key, node in aamuts['annotations'].items():
        if key == 'G': g = list(range(node['start'],node['end']))

    G_gene_length = (math.floor(len(g)/3)+1)
    N_term = [*range(1,36)]
    T_Membrane = [*range(36,67)]
    Mucinlike_I = [*range(67,164)]
    Centralconserved = [*range(164,186)] 
    Heparin_binding = [*range(186,224)]
    Mucinlike_II = [*range(224, G_gene_length+1)]
    
    dictofdomains = {'N-terminal':N_term,'Transmembrane':T_Membrane,'Mucin-like I':Mucinlike_I,'Central conserved domain':Centralconserved,'Heparin-binding domain':Heparin_binding,'Mucin-like II': Mucinlike_II}
    dictionary_synonymous =dict()
    dictionary_nonsynonymous=dict()
    for k, l in dictofdomains.items():
        list3=[]
        for i in synonymouslist:
            if i in l:list3.append(i)
        dictionary_synonymous[k] =list3
        list1=[]
        for j in a['G']:
            if j in l: list1.append(j)
        dictionary_nonsynonymous[k] =list1
    regions =[]
    nonsyn =[]
    syn =[]
    ratios=[]
    for (i, j), (k, l) in zip(dictionary_nonsynonymous.items(), dictionary_synonymous.items()):
        regions.append(i), nonsyn.append(len(j)), syn.append(len(l))
    for a, b in zip(nonsyn, syn):
        ratio = a/b
        ratios.append(ratio)
    df = pd.DataFrame({'Region': regions,'Nonsynonymous Mutations': nonsyn, 'Synonymous Mutations':syn, "dN/dS ratio": ratios })
    df.to_csv('G_gene.csv')

def non_synonymous_or_synonymous(aa_muts, nt_muts):
    aa = AA_Mutations(aa_muts, nt_muts)
    mu = MutationsineachGene(aa_muts, nt_muts)

    synonymousmutations =[]
    nonsynonymousmutations =[]
    ratios=[]
    sel =[]
    listofgenes =('F','G','M','M2','NS1','NS2','P','SH','N', 'L')
    for gene in listofgenes:
        list1 =[]
        list2 =[]
        for (k,l), (i,j) in zip(mu.items(), aa.items()):
            if k == gene and i == gene:
                a =(list(l['Codon']))
                c = Counter(a)
                b = Counter(j)
                d = c-b
                

        for i, j in b.items():
            list1.append(j)
        nonsynonymousmutations.append(sum(list1))
        for k, l in d.items():
            list2.append(l)
        synonymousmutations.append(sum(list2))
        for a, b in zip(nonsynonymousmutations, synonymousmutations):
            ratio = a/b
            if ratio>1:selection =('adaptive')
            elif ratio<1: selection = ('purifying')
            elif ratio ==1: selection =('neutral')
        sel.append(selection)
        ratios.append(ratio)

    df = pd.DataFrame({"gene":listofgenes, "synonymous mutations": synonymousmutations, "nonsynonymous mutations":nonsynonymousmutations, "dN/dS ratio":ratios })
    df.to_csv('synonymousnonsynonymous.csv')

print(non_synonymous_or_synonymous('aa_muts.json', 'nt_muts.json'))
