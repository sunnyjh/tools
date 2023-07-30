#_*_ coding: utf-8 _*
import os,re,sys
from Bio import Seq
from Bio import SeqIO
import regex as re
from Bio.SeqUtils import molecular_weight


def s2(sequence):
   ###氨基酸残基的分子量（去掉了一个水分子的质量）
   MW={
   'A':71.07,'C':103.10,'D':115.08,'E':129.11,'F':147.17,
   'G':57.05,'H':137.14,'I':113.15,'K':128.17,'L':113.15,
   'P':97.11,'N':114.10,'M':131.19,'Q':128.13,'R':156.18,
   'S':87.07,'T':101.14,'V':99.13,'W':186.20,'Y':163.17
   }
    
   sum = 0 
   for i in range(0,len(sequence)):
       if sequence[i] in MW:
           sum += (MW[sequence[i]])
   mass = (sum+18)
   return (mass)

def ORF(input_seq):
    startP = re.compile('ATG')
    nuc = input_seq.replace(' ','')
    longest = (0,)
    for m in startP.finditer(nuc, overlapped=True):
        if len(Seq.Seq(nuc)[m.start():].translate(to_stop=True)) > longest[0]:
            pro = Seq.Seq(nuc)[m.start():].translate(to_stop=True)
            longest = (len(pro),m.start(), str(pro), nuc[m.start():m.start()+len(pro)*3+3])

    return(longest)

#################################################
main_NM_file,refgene,refgeneMrna = sys.argv[1:]

main_NM={}
protein_NM={}

###
with open(main_NM_file,"r") as f:
    for line in f:
        arr = line.strip("\n").split("\t")
        if arr[1].startswith("NR"):continue 
        main_NM[arr[1]] = arr[0]

###        
gene_coord = open("gene_coord.txt","w")
with open(refgene,"r") as f:
    for line in f:
        arr = line.strip("\n").split("\t")
        if not arr[1] in main_NM: 
            continue
        if "_" in arr[2]:
            continue

        gene_coord.write("{0}\t{1}\t{2}\t{3}\n".format(arr[2],arr[4],arr[5],arr[12])) 
        protein_NM[arr[1]] = 1

gene_coord.close()
###
gene_molecular_weight = open("gene_mass.xls","w")
gene_molecular_weight.write("Gene\tTranscript\tMolecular Weight\n")
with open(refgeneMrna,"r") as f:
    for line in f:
        line = line.strip("\n")
        transcript = "" 
        if line.startswith(">"): 
            transcript = line.split(" ")[0].replace(">","")
            if transcript in protein_NM:
                mRNA_seq = next(f).strip("\n")
                #find ORF and translate
                orf = ORF(mRNA_seq)
                protein_len,protein = orf[0],orf[2]
                #other molecular mass method
                #protein_mass = s2(protein)
                protein_mass = molecular_weight(protein,"protein")
                gene_molecular_weight.write("{0}\t{1}\t{2}kDa\n".format(main_NM[transcript],transcript,str(round(protein_mass/1000,1))))

gene_molecular_weight.close()
          
        
