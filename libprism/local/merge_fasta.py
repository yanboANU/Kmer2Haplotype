#########################################################################
# File Name: merge_fasta.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Sat 18 Sep 2021 06:00:06 PM AEST
#########################################################################
#!/bin/bash
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

A_records = []
B_records = []
phased_records = []
unphased_records = []

path=sys.argv[1]
for root, dirs, files in os.walk(path):
    for d in dirs:
        #print d
        filename=path+"/"+d+"/assembly.fasta"
        if os.path.exists(filename) == False:
            print (filename, "not exist")
            continue 
        if d.startswith("groupA_") or d.startswith("groupB_"):
            for seq_re in SeqIO.parse( filename, "fasta"):
                seq_re.id =d + '_' + seq_re.id
                #print "phased",seq_re.id
                phased_records.append(seq_re)
                if d.startswith("groupA_"):
                    A_records.append(seq_re)
                if d.startswith("groupB_"):
                    B_records.append(seq_re)
        elif d.startswith("groupAB") or d.startswith("group_nonSNV"):
            for seq_re in SeqIO.parse( filename, "fasta"):
                seq_re.id =d + '_' + seq_re.id
                #print "unphased",seq_re.id
                unphased_records.append(seq_re)
          
    break


SeqIO.write(A_records, 'A_contigs.fasta', 'fasta')
SeqIO.write(B_records, 'B_contigs.fasta', 'fasta')
SeqIO.write(phased_records, 'phased_contigs.fasta', 'fasta')
SeqIO.write(unphased_records, 'unphased_contigs.fasta', 'fasta')
