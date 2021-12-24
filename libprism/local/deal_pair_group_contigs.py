#!/usr/bin/env python3.4

# goal: according blasr result
# groupA1, contig1 ----------------
# groupB1, contig1  ----- ----- contig2
# experimental show: most group pair alignment info not that clear
# when pair has clear info, it still have gaps



import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_phased_contigs(filename):
    contigsA, contigsB = {}, {}
    lA, lB =0, 0
    for seq_re in SeqIO.parse( filename, "fasta"):
        words = seq_re.id.split('_')
        #ID = words[0] + '_' + words[1]
        ID = words[1]
        if seq_re.id.startswith("groupA"):
            if ID not in contigsA:
                contigsA [ ID ]=[]
            contigsA[ID].append( seq_re )
            lA += len(seq_re.seq)
        else:
            if ID not in contigsB:
                contigsB [ ID ]=[]
            contigsB[ID].append( seq_re )
            lB += len(seq_re.seq)
      
    print ('phased A contig size:', len(contigsA) )
    print ('phased B contig size:', len(contigsB) )

    print ('phased A contig length:', lA)
    print ('phased A contig length:', lB)
    return contigsA, contigsB

if __name__ == "__main__":

    #input: phased_contigs.fasta 
    contigsA, contigsB = read_phased_contigs(sys.argv[1])

    a = contigsA.keys()-contigsB.keys()
    b = contigsB.keys()-contigsA.keys()

    print ("only in contigsA")
    Asize =0
    for ID in a:
        for seq_re in contigsA[ID]:
            temp =len(seq_re.seq)
            print (seq_re.id, temp)
            Asize += temp
    print ("A extra length", Asize) 

    print ("only in contigsB") 
    Bsize =0
    for ID in b:
        for seq_re in contigsB[ID]:
            temp =len(seq_re.seq)
            print (seq_re.id, temp)
            Bsize += temp
            
    print ("B extra length", Bsize)
    
    #os.system("mkdir group_pair_blasr") 
    for ID in contigsA:
        if ID in a:
            continue
        Asize, Bsize = 0,0   
        for seq_re in contigsA[ID]:
            temp =len(seq_re.seq)
            #print (seq_re.id, temp)
            Asize += temp
        for seq_re in contigsB[ID]:
            temp =len(seq_re.seq)
            #print (seq_re.id, temp)
            Bsize += temp
        if (abs(Asize-Bsize)>0.1*min(Asize, Bsize)):  
            print (ID, Asize, Bsize)
        '''
        if len(contigsA[ID]) != len(contigsB[ID]):
            print (ID, "two groups diff")
            for seq_re in contigsA[ID]:
                print (seq_re.id, len(seq_re.seq))
            for seq_re in contigsB[ID]:
                print (seq_re.id, len(seq_re.seq))
               
            command = "/g/data/te53/software/MaSuRCA/MaSuRCA-3.4.2/CA8/Linux-amd64/bin/blasr "   
            command = command + "groupA_" + ID + "/assembly.fasta" + " groupB_" + ID 
            command = command + "/assembly.fasta" + " -nproc 20 -m 4 >group_pair_blasr/" +ID +".blasr" 
            print (command)
            os.system(command)
            #sys.exit()
        '''
           
