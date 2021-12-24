#!/usr/bin/env python3.4

import os
import sys
import tools

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_blasr(filename):
    ######
    # if contigs align two sym refs, keep both alignment
    # otherwise, only keep best alignment
    #####
    contigMap = {}
    refMap = {}
    dealedC =set()
    partAlign = set() # contigID set
    with open(filename, "r") as f:
        for line in f:
            words = line.strip().split(' ')
            contigID=words[0]
            if contigID in dealedC:
                continue
            cLen = int(words[1])
            if cLen < 5000:
                continue
            cStart = int (words[2])
            cEnd = int (words[3])
            refID  = words[6]
            rLen = int(words[7])
            rStart = int (words[8])
            rEnd = int (words[9])
            numMatch = int(words[12]) 
            
            if contigID == refID:
                continue
            if numMatch < 2000 or numMatch < 0.1*max(cLen, rLen) :
                continue
            if rStart> 0.01*rLen and rEnd< 0.99*rLen and cStart> 0.01*cLen and cEnd< 0.99*cLen:
                partAlign.add(contigID)
                continue

            score = int(words[11]) 
            if contigID not in contigMap:
                contigMap[ contigID ] = []
                contigMap[ contigID ].append(  words )
    return contigMap, refMap, partAlign

def read_contigs(filename):
    contigs = {} 
    length = 0
    shortC = set()
    for seq_re in SeqIO.parse( filename, "fasta"):
        contigs[ seq_re.id ] = seq_re
        l = len(seq_re.seq)
        if l <5000:
            shortC.add(seq_re.id)
        length += l
    return contigs, shortC, length


def deal_covered_unphased_contigs(contigMap):
    ######### 
    # remove contigs covered by phased refs
    ##########
    inSize, inNum = 0, 0
    dealed = set() # dealed covered unphased contigID
    removeSet= set()
    for c in contigMap:
        l = len( contigMap[c] )
        #print (r, l)
        for i in range(l):
            words = contigMap[c][i]
            cID = words[0] 
            cLen = int(words[1])
            cStart = int (words[2])
            cEnd = int (words[3])
            cD = words[4] # direction

            rName = words[6]
            rLen = int (words[7])
            rStart = int (words[8])
            rEnd = int (words[9])
            rD = words[10] # direction

            numMatch = int(words[12])
            #if cStart <= 100 and cEnd >= cLen-100: #ignore extend length smaller than 200
            if (cStart <= 0.1*cLen and cEnd >= cLen*0.9):
                print ("covered")
                print (words[:11])
                if cID not in dealed:
                    inSize += cLen
                    inNum += 1
                    #print (cID, "alreay covered by other phased contigs")
                    dealed.add(cID)
                removeSet.add( cID )
            elif (cD==rD and rStart>=cStart and rLen-rEnd>=cLen-cEnd):

                overhang = cStart
                overhang += cLen - cEnd
                if overhang > numMatch:
                    continue
                print ("++ covered")
                print (words[:11])
                if cID not in dealed:
                    inSize += cLen
                    inNum += 1
                    #print (cID, "alreay covered by other phased contigs")
                    dealed.add(cID)
                removeSet.add( cID )
            elif (cD!=rD and rStart>=cLen-cEnd and cStart<=rLen-rEnd):

                overhang = cStart
                overhang += cLen - cEnd
                if overhang > numMatch:
                    continue
                print ("+- covered")
                print (words[:11])
                if cID not in dealed:
                    inSize += cLen
                    inNum += 1
                    #print (cID, "alreay covered by other phased contigs")
                    dealed.add(cID)
                removeSet.add( cID )

    print (inNum, "unknow contigs are covered by other unknow contigs")
    print (inSize, "length unphased contigs are covered by phased contigs")
    return dealed



if __name__ == "__main__":

    #input: blasr, unknow.fasta
    filename = sys.argv[1]
    contigMap, refMap, partAlign = read_blasr(filename)
    un_contigs, un_s, l1 = read_contigs(sys.argv[2]) # un_s: unphased short contigs
    print ('unkonw contig size:', len(un_contigs) )
    print ('unknow contig length:', l1)


    coveredSet = deal_covered_unphased_contigs(contigMap) #unphased contigs are covered by phased contig
