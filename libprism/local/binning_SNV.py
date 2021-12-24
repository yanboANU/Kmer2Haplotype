#########################################################################
# File Name: binning_TGS_reads.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Fri 07 May 2021 14:53:44 AEDT
#binning snv by voting 
#########################################################################
#!/bin/bash
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tools

#if all reads are phased corrected,
#according to phased reads, get phased SNV
#import build_matrix_for_supperReads


def write_phased_ID(A, B, filename):

    with open(filename, "w") as fout:
        fout.write("group_1_A_groupsize_%s\n" % len(A)) 
        for ele in A:
            fout.write("%s " % ele)
        fout.write("\n") 
        fout.write("group_1_B_groupsize_%s\n" % len(B)) 
        for ele in B:
            fout.write("%s " % ele)
        fout.write("\n") 
    return 

def read_matrix(filename, group1, group2):

    node2Group1 = {}
    node2Group2 = {}
    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            wordsLen = len( words )
            nodeNumber = int( words[0] )
            readID = int(words[1].split('S')[1].split('_')[0])
            #print (readID)
            #sys.exit()

            if nodeNumber < 0:
                continue
            
            if readID % 2 == 1:
                for i in range(2, wordsLen-1, 2):
                    nodeID = words[i]
                    if nodeID not in node2Group1:
                        node2Group1[ nodeID ] = 0
                    node2Group1[ nodeID ] += 1
            else:
                for i in range(2, wordsLen-1, 2):
                    nodeID = words[i]
                    if nodeID not in node2Group2:
                        node2Group2[ nodeID ] = 0
                    node2Group2[ nodeID ] += 1
     
    print (len(node2Group1) )
    print (len(node2Group2) )
    unphasedNode =set()
    for nodeID in node2Group1:
        if nodeID in node2Group2:
            if node2Group1[ nodeID ] > node2Group2[ nodeID]:
                group1.add(nodeID)
            elif node2Group1[ nodeID ] < node2Group2[ nodeID] :
                group2.add(nodeID)
            else:
                print (nodeID, node2Group1[ nodeID ], node2Group2[ nodeID])
                unphasedNode.add(nodeID)
        else:
            group1.add(nodeID)

    for nodeID in node2Group2: 
        if (nodeID not in node2Group1):
            group2.add(nodeID)

    print ("groupA size", len(group1))
    print ("groupB size", len(group2))
    print ("unphased number", len(unphasedNode))
    for n in unphasedNode:
        symN = tools.get_symmetrical_node(n)
        if symN in group1:
            group2.add(n)
        if symN in group2:
            group1.add(n)
    print ("groupA size", len(group1))
    print ("groupB size", len(group2))
    return group1, group2





# first can try build_matrix_for_supperReads.py
if __name__ == '__main__':
    
    if len(sys.argv) < 2:
        print ("input: TGS_matrix, output: phased_kmer_ID")
        sys.exit()

    groups1, groups2 = set(), set()

    print ("reads TGS matrix")
    A, B = read_matrix(sys.argv[1], groups1, groups2)
    write_phased_ID(A, B, "phased_kmer_ID")
