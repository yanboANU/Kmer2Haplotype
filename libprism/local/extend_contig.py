#########################################################################
# File Name: extend_contig.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Mon 31 May 2021 05:42:34 PM AEST
#########################################################################
#!/bin/bash
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_fasta(filename):
    records = {}
    #recordIDs = set()
    for seq_re in SeqIO.parse(filename, "fasta"):
        #print(seq_re.id)
        records [seq_re.id] = seq_re 
        #recordIDs.add( seq_re.id )
    return records #, recordIDs 

def write_fasta(records):
    os.system('pwd')
    print ("write merged fasta")
    recordList = []
    for ID in records:
        recordList.append(records[ID])
    SeqIO.write(recordList, "merged.fasta", "fasta")
    return


def read_blasr(filename):

    nextC, preC = {}, {}
    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            if len(words) == 19:
                refName = words[0]
                refLen = int(words[1])
                start, end = int(words[2]), int(words[3])
                assert start < end
                queryName = words[5]

                if ( (refName.startswith('groupA') and queryName.startswith('groupB')) or
                     (refName.startswith('groupB') and queryName.startswith('groupA') ) ):
                    refGroupID = refName.split('_')[1]
                    qGroupID = queryName.split('_')[1]
                    if refGroupID == qGroupID:
                        #print (refName, queryName, "groupID", refGroupID, qGroupID)
                        continue
                refSeq = words[-3]
                matchSeq = words[-2]
                querySeq = words[-1]
                # contig must align 80\%
                queryLen = int( words[6] )
                queryAlignStart, queryAlignEnd = int(words[7]), int(words[8])
                assert queryAlignStart < queryAlignEnd
                if refName == queryName:
                    continue
                if (end-start)/float(refLen) > 0.9 or end-start<2000:
                    continue
                if (queryAlignEnd-queryAlignStart)/float(queryLen) > 0.9:
                    continue
                if ( (start < 1000 or queryAlignStart < 1000) 
                        and (end > refLen - 1000 or queryAlignEnd > queryLen - 1000) ):
                    
                    numMatch, numMismatch = int(words[11]), int(words[12])
                    numIns, numDel = int(words[13]), int(words[14])
                    print (refName, refLen, start, end, queryName, queryLen, queryAlignStart, queryAlignEnd)
                    print (numMatch, numMismatch, numIns, numDel)
                    rate =  (numMismatch + numIns + numDel )/float(numMatch) *100
                    print (rate)
                    if rate > 0.2:
                        continue
                    if start < 1000:
                        
                        if refName not in preC:
                            preC[ refName ] = [] #set()
                        preC[refName].append( words[:16]  )    
                        ''' 
                        if queryName not in nextC:
                            nextC[ queryName ] = [] #set()
                        nextC[queryName].append( words[:16] )    
                        '''
                    elif queryAlignStart < 1000:
                        if refName not in nextC:
                            nextC[ refName ] = [] #set()
                        nextC[refName].append( words[:16] )
                        '''  
                        if queryName not in preC:
                            preC[ queryName ] = [] #set()
                        preC[queryName].append( words[:16]  )    
                        '''
    for c in preC:
        print (c, preC[c])
    for c in nextC:
        print (c, nextC[c])
      
    return preC, nextC

def extend_seqs(preC, nextC, seqs):

    print (type(seqs))
    seqsKey = seqs.keys()
    removeR = set() # remove records, after merge
    for c in nextC:
        print (c, nextC[c])
        if len(nextC[c]) == 1:
            #'groupA_41_contig_1', '313630', '305144', '313630', '+', 
            #'group_nonSNV_contig_31', '18726', '0', '8488', '+', '-42342', '8480', '3', '3', '5', '0'
            #qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand 
            #score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq 
            eles = nextC[c][0]
            c1, l1, s1, e1, d1, c2, l2, s2, e2, d2 = eles[:10]
            score, numMatch, numMismatch, numIns, numDel, mapQV = eles[10:]
            print ( c1, l1, s1, e1, d1, c2, l2, s2, e2, d2 )
            print (score, numMatch, numMismatch, numIns, numDel, mapQV)
            assert c == c1
            if c2 in preC: 
                if len(preC[c2]) == 1:
                    assert c1 in seqs 
                    assert c2 in seqs
                    print (seqs[c1].id, seqs[c2].id)
                    #sys.exit()
                    if d1 == d2:
                        print seqs[c1].seq[int(s1): int(e1)]
                        print seqs[c2].seq[int(s2): int(e2)]
                        record = SeqRecord(seqs[c1].seq + seqs[c2].seq[int(e2):] , id=seqs[c1].id+"_"+ seqs[c2].id )
                        print (record.id)
                        print (len(record.seq))
                        removeR.add( seqs[c1].id )
                        removeR.add( seqs[c2].id )
                        seqs[ record.id ] = record
                        #sys.exit()
                    else:
                        print ("align different direction") 
                        print (seqs[c1].id, seqs[c2].id)
                        print seqs[c1].seq[int(s1): int(e1)]
                        print seqs[c2].seq[int(s2): int(e2)]
        else:
            print ("")

    for ID in removeR:
        seqs.pop(ID)
    return


#input: *.blasr *fasta

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print ("*blasr *fasta")
        sys.exit()
    preC, nextC = read_blasr(sys.argv[1])
    seqs = read_fasta(sys.argv[2])
    extend_seqs(preC, nextC, seqs)
    write_fasta(seqs) # write final results
