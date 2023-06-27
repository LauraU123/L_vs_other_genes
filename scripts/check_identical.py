from Bio import SeqIO
import argparse

def check_identical(linput, restinput, loutput, restoutput):
    L_input = SeqIO.parse(linput, "fasta")
    rest_input = SeqIO.parse(restinput, "fasta")
    L_entries = []
    rest_entries= []
    L_output = []
    rest_output = []
    for entry in L_input: 
        L_entries.append(entry.id)
    for entry in rest_input: 
        rest_entries.append(entry.id)
    common = [value for value in L_entries if value in rest_entries]  
    print(common)


    L_input = SeqIO.parse(linput, "fasta")
    rest_input = SeqIO.parse(restinput, "fasta")
    for entry in L_input: 
        print(entry.id)
        if entry.id in common:  
            print(entry)
            L_output.append(entry)
            
    for entry in rest_input: 
        print(entry.id)
        if entry.id in common:  
            rest_output.append(entry)

    SeqIO.write(rest_output, restoutput, "fasta")
    SeqIO.write(L_output, loutput, "fasta")        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="make new reference depending on whether the entire genome or only part is to be used for the tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--inputl", required=True, help="fasta input with alignment of l")
    parser.add_argument("--inputrest", required=True, help="fasta input with alignment of rest of the genome")
    parser.add_argument("--outputl", required=True, help="fasta output with alignment of l")
    parser.add_argument("--outputrest", required=True, help="fasta output with alignment of rest")
    args = parser.parse_args()

    check_identical(args.inputl, args.inputrest, args.outputl, args.outputrest)