#########################################################################
# File Name: using_extend_kmer.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Wed 04 Dec 2019 11:44:21 AEDT
#########################################################################
#!/bin/bash
import os
import sys
#import  tools
from libprism.local import tools
#import build_matrix_for_PairEndReads

def assign_ID_2_kmerPair(kmerPair):
    
    heteKmer = []
    heteID = []
    count = 1
    for (k1, k2, w) in kmerPair:
        ID1 = str(count) + '_0'
        ID2 = str(count) + '_1'
        count += 1
        heteKmer.append( (k1, ID1) )
        heteKmer.append( (k2, ID2) )
        heteID.append( (ID1, ID2) )
    return heteID, heteKmer

def read_kmer_pair(filename):
    
    kmerPair = []
    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            kmerPair.append(words)
    return kmerPair

def scan_seq(seq, k):

    kmers = []
    seqLen = len(seq)
    #for i in range(seqLen - k):
    for i in range(seqLen - k + 1): #debug
        key = str(seq[i:i + k])
        if key.count('A') + key.count('C') + key.count('T') + key.count('G')  != k:
            continue
        Rkey = tools.reverse(key)
        if key < Rkey:
            kmers.append(key)
        else:
            kmers.append(Rkey)

    return kmers   

def get_ambigous_kmer_ID(extendKmers):
    
    ambigousID=set()  
    for kmer in extendKmers:
        if len(extendKmers[kmer]) > 1:
            #print("%s %s" % (kmer, extendKmers[kmer]) )
            for ele in extendKmers[kmer]:
                ambigousID.add(ele)
                ambigousID.add( tools.get_symmetrical_node(ele) )

    print ("ambigousID number", len(ambigousID) )
    print ("ambigousID", ambigousID)
    return ambigousID

def assign_ID2_extendKmer(extendSNPPair, heteID, k):

    heteIDLen = len(heteID)
    print ("kmer size", k)
    assert heteIDLen == len(extendSNPPair)
    
    extendKmers= {}
    for i in range(heteIDLen):
        ek1, ek2, w = extendSNPPair[i]
        ID1, ID2 = heteID[i]
        kmers1 = scan_seq(ek1, k)
        kmers2 = scan_seq(ek2, k)
        for kmer in kmers1:
            if kmer not in extendKmers:
                extendKmers [ kmer] = []
            extendKmers [ kmer ].append(ID1)
            #fout.write("%s %s\n" % (kmer, ID1) )
        
        for kmer in kmers2:
            if kmer not in extendKmers:
                extendKmers [ kmer] = []
            extendKmers [ kmer ].append(ID2)
            #fout.write("%s %s\n" % (kmer, ID2) )
    return extendKmers

def write_heterozygous_snp_kmer(heteKmer, extendKmers, ambigousID):

    with open("kmer_ID_mid", "w") as fout:
        for (kmer, ID) in heteKmer:
            if ID not in ambigousID:
                fout.write("%s %s\n" % (kmer, ID) )
           
    print ("after assign extend kmer, one kmer multiple IDs")       
    with open("kmer_ID", "w") as fout:
        for kmer in extendKmers:
            if len(extendKmers[kmer]) == 1 and extendKmers[kmer][0] not in ambigousID:
                fout.write("%s %s\n" % (kmer, extendKmers[kmer][0]) )
            #else:
                #print("%s %s" % (kmer, extendKmers[kmer]) )
    return

def get_heterozygous_snp_kmer(file1, file2):

    SNPPair = read_kmer_pair( file1 ) 
    extendSNPPair = read_kmer_pair( file2 )
    heteID, heteKmer = assign_ID_2_kmerPair(SNPPair)
    k = len(SNPPair[0][0])
    extendKmers = assign_ID2_extendKmer(extendSNPPair, heteID, k)
    ambigousID = get_ambigous_kmer_ID(extendKmers)
    write_heterozygous_snp_kmer(heteKmer, extendKmers, ambigousID)

    return

'''
if __name__ == '__main__':
    
    if len(sys.argv) < 3:

        print ("python k_21_pair.snp k_21_pair.snp.extend (must same order)")
        print ("python k_21_pair.snp k_21_pair.snp.extend kmer_pair_file")
        print ("python k_21_pair.snp k_21_pair.snp.extend h1_21mers h2_21mers")
        sys.exit()
        
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    SNPPair = read_kmer_pair( file1 ) 
    extendSNPPair = read_kmer_pair( file2 )
    print ("paramerter number", len(sys.argv)-1)
    if len(sys.argv) - 1 == 2:
        heteID, heteKmer = assign_ID_2_kmerPair(SNPPair)
    elif len(sys.argv) - 1 == 4:
        heteID, heteKmer =  assign_ID_2_kmerPair_based_ground_truth(SNPPair, sys.argv[3], sys.argv[4])
    elif len(sys.argv) - 1 == 3:   #/home/yanbo/bio/Kmer2Haplotype/realData/A.thaliana/FUNZIP/cvi0-col0-SNPs/cvi_col_kmer_pair_k31
        heteID, heteKmer =  assign_ID_2_kmerPair_based_one_file(SNPPair, "/home/yanbo/bio/Kmer2Haplotype/realData/A.thaliana/FUNZIP/cvi0-col0-SNPs/cvi_col_kmer_pair_k31")

    k = len(SNPPair[0][0])
    extendKmers = assign_ID2_extendKmer(extendSNPPair, heteID, k)
    ambigousID = get_ambigous_kmer_ID(extendKmers)
    write_heterozygous_snp_kmer(heteKmer, extendKmers, ambigousID)

def assign_ID_2_kmerPair_based_ground_truth(kmerPair, h1_file, h2_file):
    
    file1 = "/home/yanbo/bio/Kmer2SNP/experiment/realData/fungi_c.albicans/two-haplotypes/alignment_and_pair_kmer/" + h1_file
    file2 = "/home/yanbo/bio/Kmer2SNP/experiment/realData/fungi_c.albicans/two-haplotypes/alignment_and_pair_kmer/" + h2_file
    kmer2ID_0 = build_matrix_for_PairEndReads.read_kmer_ID(file1)
    kmer2ID_1 = build_matrix_for_PairEndReads.read_kmer_ID(file2)
    heteKmer = []
    heteID = []

    count = 1
    for (k1, k2, w) in kmerPair:
        flag = False
        if (k1 in kmer2ID_0) and (k2 in kmer2ID_1):
            ID1 = kmer2ID_0[k1]
            ID2 = kmer2ID_1[k2]
            if ID1[:-1] == ID2[:-1]: # k1 and k2 may from different chromosome 
                flag = True
        elif (k1 in kmer2ID_1) and (k2 in kmer2ID_0):
            ID1 = kmer2ID_1[k1]
            ID2 = kmer2ID_0[k2]
            if ID1[:-1] == ID2[:-1]: 
                flag = True
        if flag == False:
            ID1 = str(count) + '_0'
            ID2 = str(count) + '_1'
            count += 1
        heteKmer.append( (k1, ID1) )
        heteKmer.append( (k2, ID2) )
        heteID.append( (ID1, ID2) )
    return heteID, heteKmer

def assign_ID_2_kmerPair_based_one_file(kmerPair, filename):
    
    kmer2ID = read_pair_kmer_ID(filename)
    heteKmer = []
    heteID = []

    count = 1
    for (k1, k2, w) in kmerPair:
        if ( (k1 in kmer2ID) and (k2 in kmer2ID) and 
        len(kmer2ID[k1]) == 1 and len(kmer2ID[k2]) == 1 ):
            ID1 = kmer2ID[k1][0]
            ID2 = kmer2ID[k2][0]
            #print (ID1, ID2)
            #sys.exit()
            #if ID1[:-1] == ID2[:-1]: # k1 and k2 may from different chromosome  
        else:
            ID1 = str(count) + '_0'
            ID2 = str(count) + '_1'
            count += 1
        heteKmer.append( (k1, ID1) )
        heteKmer.append( (k2, ID2) )
        heteID.append( (ID1, ID2) )
    return heteID, heteKmer

def read_pair_kmer_ID(filename): # ground truth

    kmer2ID = {} 
    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            k1, k2 = words[0], words[1]
            h1, h2 = words[-2], words[-1]
            if k1 not in kmer2ID:
                kmer2ID [k1] = []

            if k2 not in kmer2ID:
                kmer2ID[ k2 ] = []
            kmer2ID [ k1 ].append( words[2] + "_" + h1 )
            kmer2ID [ k2 ].append( words[2] + "_" + h2 )
 
    print ("one kmer multiple ID")
    count = 0
    for key in kmer2ID:
        if len(kmer2ID[key]) > 1:
            #print (key, kmer2ID[key])
            count += len(kmer2ID[key])
    print ("end")
    print ("%s kmers not unique kmer on two haplotypes" % count)
    return kmer2ID

'''
