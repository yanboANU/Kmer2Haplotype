########################################################################
# File Name: binning_TGS_reads.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Tue 03 Dec 2019 14:53:44 AEDT
#binning TGS reads by voting 
# last version in ../local-1
#########################################################################
#!/bin/bash
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tools

#import build_matrix_for_supperReads


def read_phased_ID(filename):

    groups1 = []
    groups2 = []
    count = 1
    with open(filename, "r") as f:
        for line in f:
            if count%4 == 2:
                words = line.split()
                if len(words)<=30: # remove small phasedID group
                    continue
                else:
                    groups1.append(words)

            if count%4 == 0:
                words = line.split()
                if len(words)<=30:
                    continue
                else:
                    groups2.append(words)
            count += 1
    print ("groups1 size", len(groups1))       
    print ("groups2 size", len(groups2))      
    for g in groups1:
        print (len(g), end=' ')
    print ('')    
    for g in groups2:
        print (len(g), end=' ')

    print ('')    
    return groups1, groups2

def merge_group(readID, inter1, inter2, m):
    
    nonZeroIndex1 = tools.get_non_zero_index(inter1)
    nonZeroIndex2 = tools.get_non_zero_index(inter2)
    if len(nonZeroIndex1) + len(nonZeroIndex2) == 2:

        if len(nonZeroIndex1) == 2:
            pos1 = nonZeroIndex1[0]
            pos2 = nonZeroIndex1[1]

            n1 = str(pos1) + "_1"
            n2 = str(pos2) + "_1"
            #print ("merge ", n1, n2, inter1[pos1], inter1[pos2])

            if (n1, n2) not in m:
                m[(n1, n2)] = []
            m[(n1, n2)].append( (inter1[pos1], inter1[pos2] ) )    

        elif len(nonZeroIndex2) == 2:
            pos1 = nonZeroIndex2[0]
            pos2 = nonZeroIndex2[1]

            n1 = str(pos1) + "_0"
            n2 = str(pos2) + "_0"
            #print ("merge ", n1, n2, inter2[pos1], inter2[pos2])

            if (n1, n2) not in m:
                m[(n1, n2)] = []
            m[(n1, n2)].append( (inter2[pos1], inter2[pos2] ) )   

        if len(nonZeroIndex1) == 1 and len(nonZeroIndex2) == 1:
            pos1 = nonZeroIndex1[0]
            pos2 = nonZeroIndex2[0]
            if pos1 != pos2: #and inter1[pos1] == inter2[pos2]:

                n1 = str(pos1) + "_1"
                n2 = str(pos2) + "_0"
                #print ("merge ", n1, n2, inter1[pos1], inter2[pos2])

                if (n1, n2) not in m:
                    m[(n1, n2)] = []
                m[(n1, n2)].append( (inter1[pos1], inter2[pos2] ) )   
    return


def update_read_set(readID, inter1, inter2, A, B, C):
   
    #############
    # only inter with A[i] or B[i] --> True; inter with one A[i] and one B[i] and A[i]><B[i] -->True
    # otherwise False
    ############
    nonZeroIndex1 = tools.get_non_zero_index(inter1)
    nonZeroIndex2 = tools.get_non_zero_index(inter2)
    if len(nonZeroIndex1) <= 1 and len(nonZeroIndex2) <= 1: 
        if len(nonZeroIndex1) == 1 and len(nonZeroIndex2) == 0:
            pos = nonZeroIndex1[0]
            if pos not in A:
                A[pos] = set()
                B[pos] = set()
            A[pos].add(readID)     
        elif len(nonZeroIndex1) == 0 and len(nonZeroIndex2) == 1:
            pos = nonZeroIndex2[0]
            if pos not in A:
                A[pos] = set()
                B[pos] = set()
            B[pos].add(readID)   
        elif len(nonZeroIndex1) == 1 and len(nonZeroIndex2) == 1:
            pos1 = nonZeroIndex1[0]
            pos2 = nonZeroIndex2[0]
            if pos1 == pos2:
                if inter1[pos1] > inter2[pos2]:    
                    if pos1 not in A:
                        A[ pos1 ] = set()
                        B[ pos1 ] = set()
                    A[ pos1 ].add(readID)  
                elif inter1[pos1] < inter2[pos2]:    
                    if pos2 not in B:
                        A[ pos2 ] = set()
                        B[ pos2 ] = set()
                    B[ pos2 ].add(readID)  
                else:
                    C.add(readID)
                    return False
            else:
                C.add(readID)
                return False
    else:            
        C.add(readID)
        return False

    return True
def read_matrix(filename, groups1, groups2, m, k):

    size = len(groups1)
    A, B = {}, {}
    C = set() # ambiguous readID

    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            wordsLen = len( words )
            nodeNumber = int( words[0] )
            readID = words[1]
            if nodeNumber < 0:
                continue
            nodeID = []
           
            #print (line)
            for i in range(2, wordsLen-1, 2):
                supportNum = int(words[i+1].split('_')[0])
                #print (supportNum)
                if k<31 and supportNum <= 1:
                    continue
                nodeID.append(words[i]) # with ground truth
            #sys.exit()  
            
            if k<31 and len(nodeID) <=1:
                continue
               
            #try_merge_block_according_TGS_reads(readID, nodeID, groups1, groups2)
            
            inter1, inter2 = [], []
            for i in range(size):
                l1 = set( nodeID ).intersection( set(groups1[i]) ) 
                inter1.append( len(l1) )
            
            for i in range(size):
                l2 = set( nodeID ).intersection( set(groups2[i]) ) 
                inter2.append( len(l2) ) 
              
            label = update_read_set(readID, inter1, inter2, A, B, C)
            
            #can be deleted, only show read inter with multiple group
            '''
            if label == False:
                merge_group(readID, inter1, inter2, m)
            '''    

    print ("phased group size", len(A), len(B))
    for i in sorted(A.keys()):        
        print ("debug A number", i, len(A[i]))
        print ("debug B number", i, len(B[i]))
   
    '''
    os.system("mkdir readsPhasing")
    foutAB = open("readsPhasing/ambiguous_readID","w")
    for readID in C:
        foutAB.write( "%s\n" % (readID) )
    foutAB.close()
    '''
    return A, B, C


def write_two_groups_reads(A, B, filenames):
    foutA = open("groupA.fastq", "w")
    foutB = open("groupB.fastq", "w")
    foutAID= open("readsPhasing/groupA_readsID", "w")
    foutBID= open("readsPhasing/groupB_readsID", "w")

    # reads not assign to any group, then put all 
    cnt, cntA, cntB = 0, 0, 0
    for f in filenames:
        print ("check", f)
        for seqRecord in SeqIO.parse(f, "fastq"):
            if seqRecord.id in A:
                cntA += 1
                tools.write_one_fastq_record(foutA, seqRecord)
                foutAID.write("%s\n" % seqRecord.id )
            elif seqRecord.id in B:
                cntB += 1
                tools.write_one_fastq_record(foutB, seqRecord)
                foutBID.write("%s\n" % seqRecord.id )
            else:
                cnt += 1
                tools.write_one_fastq_record(foutA, seqRecord)
                tools.write_one_fastq_record(foutB, seqRecord)
    print ("groupA contain %s reads, phased reads %s" % (cntA+cnt, cntA))
    print ("groupB contain %s reads, phased reads %s" % (cntB+cnt, cntB))
    print ("unphased reads %s" % (cnt))
    foutA.close()
    foutB.close()
    foutAID.close()
    foutBID.close()
    

def write_two_groups_fasta(A, B, filenames):
    foutA = open("groupA.fasta", "w")
    foutB = open("groupB.fasta", "w")
    foutAID= open("readsPhasing/groupA_readsID", "w")
    foutBID= open("readsPhasing/groupB_readsID", "w")

    # reads not assign to any group, then put all 
    cnt, cntA, cntB = 0, 0, 0
    for f in filenames:
        print ("check", f)
        for seqRecord in SeqIO.parse(f, "fasta"):
            if seqRecord.id in A:
                cntA += 1
                tools.write_one_fasta_record(foutA, seqRecord)
                foutAID.write("%s\n" % seqRecord.id )
            elif seqRecord.id in B:
                cntB += 1
                tools.write_one_fasta_record(foutB, seqRecord)
                foutBID.write("%s\n" % seqRecord.id )
            else:
                cnt += 1
                tools.write_one_fasta_record(foutA, seqRecord)
                tools.write_one_fasta_record(foutB, seqRecord)
    print ("groupA contain %s reads, phased reads %s" % (cntA+cnt, cntA))
    print ("groupB contain %s reads, phased reads %s" % (cntB+cnt, cntB))
    print ("unphased reads %s" % (cnt))
    foutA.close()
    foutB.close()
    foutAID.close()
    foutBID.close()
    

def write_multiple_groups_reads(A, B, filenames):

    foutA, foutB = {}, {}
    numA, numB = {}, {} 
    os.system("mkdir group_fastq")
    
    assert len(A) == len(B) 
    for i in A:
        fout = open("group_fastq/groupA_" + str(i) +".fastq", "w")
        foutA[i] = fout
        fout = open("group_fastq/groupB_" + str(i) +".fastq", "w")
        foutB[i] = fout
    #foutAID= open("readsPhasing/groupA_readsID", "w")
    #foutBID= open("readsPhasing/groupB_readsID", "w")

    assert len(foutA) == len(foutB)
    assert len(A) == len(foutA)
    cnt, cntA, cntB = 0, 0, 0
    readID = set()
    for f in filenames:
        print ("input fastq", f)
        for seqRecord in SeqIO.parse(f, "fastq"):
            readID.add(seqRecord.id)
            for i in A:
                if i not in numA:
                    numA[i] = 0
                    numB[i] = 0
                if seqRecord.id in A[i]:
                    cntA += 1
                    numA[i] += 1
                    tools.write_one_fastq_record(foutA[i], seqRecord)
                    break
                elif seqRecord.id in B[i]:
                    cntB += 1
                    numB[i] += 1
                    tools.write_one_fastq_record(foutB[i], seqRecord)
                    break
    print ("groupA contain %s reads, phased reads %s" % (cntA, cntA))
    print ("groupB contain %s reads, phased reads %s" % (cntB, cntB))

    for i in foutA:
        print ("%s group, groupA reads number %s, groupB reads number %s" % (i, numA[i], numB[i]) )
        foutA[i].close()
    for i in foutB:    
        foutB[i].close()
    return readID 

def write_multiple_groups_fasta(A, B, filenames):

    foutA, foutB = {}, {} 
    os.system("mkdir group_fasta")
    for i in A:
        fout = open("group_fasta/groupA_" + str(i) +".fasta", "w")
        foutA[i] = fout

    for i in B:
        fout = open("group_fasta/groupB_" + str(i) +".fasta", "w")
        foutB[i] = fout
    cnt, cntA, cntB = 0, 0, 0
    readID = set()
    for f in filenames:
        print ("input fasta", f)
        for seqRecord in SeqIO.parse(f, "fasta"):
            readID.add(seqRecord.id)
            for i in A:
                if seqRecord.id in A[i]:
                    cntA += 1
                    tools.write_one_fasta_record(foutA[i], seqRecord)
                    break
                elif seqRecord.id in B[i]:
                    cntB += 1
                    tools.write_one_fasta_record(foutB[i], seqRecord)
                    break
    print ("groupA contain %s reads, phased reads %s" % (cntA, cntA))
    print ("groupB contain %s reads, phased reads %s" % (cntB, cntB))

    for f in foutA:
        f.close()
    for f in foutB:    
        f.close()
    
    return readID   



def write_ambiguous_fasta(C, filenames):

    readsID = set()
    fout = open("group_fasta/groupAB.fasta", "w") 
    cnt = 0
    for f in filenames:
        print ("input fasta", f)
        for seqRecord in SeqIO.parse(f, "fasta"):
            readsID.add(seqRecord.id)
            if seqRecord.id in C:
                cnt += 1
                tools.write_one_fasta_record(fout, seqRecord)
    print ("ambiguous reads number %s" % (cnt))

    fout.close() 
    return readsID     

def write_ambiguous_reads(C, filenames):
 
    readsID = set()
    fout = open("group_fastq/groupAB.fastq", "w") 
    cnt = 0
    for f in filenames:
        print ("input fastq", f)
        for seqRecord in SeqIO.parse(f, "fastq"):
            readsID.add(seqRecord.id)
            if seqRecord.id in C:
                cnt += 1
                tools.write_one_fastq_record(fout, seqRecord)
    print ("ambiguous reads number %s" % (cnt))

    fout.close() 
    return readsID   


def write_nonSNV_fasta(C, filenames):

    fout = open("group_fasta/group_nonSNV.fasta", "w") 
    cnt = 0
    for f in filenames:
        print ("input fasta", f)
        for seqRecord in SeqIO.parse(f, "fasta"):
            if seqRecord.id in C:
                cnt += 1
                tools.write_one_fasta_record(fout, seqRecord)
    print ("non SNV reads number %s" % (cnt))
    fout.close() 
    return  

def write_nonSNV_reads(C, filenames):
 
    fout = open("group_fastq/group_nonSNV.fastq", "w") 
    cnt = 0
    for f in filenames:
        print ("input fastq", f)
        for seqRecord in SeqIO.parse(f, "fastq"):
            if seqRecord.id in C:
                cnt += 1
                tools.write_one_fastq_record(fout, seqRecord)
    print ("non SNV reads number %s" % (cnt))
    fout.close() 
    return 

def remove_small_group(A, B):
    
    for i in sorted(A.keys()):        
        lA = len(A[i]) 
        lB = len(B[i])
        print ("debug A number", i, lA)
        print ("debug B number", i, lB)
        if lA < 500 or lB < 500: # 20k*500=10M
            A.pop(i)
            B.pop(i)
   
    print ("after remove small size")
    print ("phased group size", len(A), len(B))
    return      

def run(phasedKmerFile, matrixFile, TGSFile, k):

    print ("reads groups")
    groups1, groups2 = read_phased_ID(phasedKmerFile) # groups contain each phasedID group

    # wait to do 
    # remove small phasedID group

    print ("reads TGS matrix")
    m={}
    A, B, C = read_matrix(matrixFile, groups1, groups2, m, k) # A contain phased reads group
  
    #for key in m:
    #    print (key, m[key])
    remove_small_group(A,B) # remove small reads group

    print ("write fastq")
    filenames = TGSFile.split(',') # multiple TGS reads file
    print ("check", filenames)

    if filenames[0].endswith("fastq"):
        allReadsID = write_multiple_groups_reads(A, B, filenames)
        D = allReadsID 
        for i in A:
            D = D - A[i]
            D = D - B[i]
        allReadsID = write_ambiguous_reads(D, filenames)
    elif filenames[0].endswith("fasta"):
        allReadsID = write_multiple_groups_fasta(A, B, filenames) 
        D = allReadsID 
        for i in A:
            D = D - A[i]
            D = D - B[i] 
        allReadsID = write_ambiguous_fasta(D, filenames)

# first can try build_matrix_for_supperReads.py
if __name__ == '__main__':
    
    if len(sys.argv) < 4:
        print ("phased_kmer_ID, TGS_matrix, sd.fastq, k")
        sys.exit()

    print (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    run(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))

    '''
    if filenames[0].endswith("fastq"):
        ####write_two_groups_reads(A, B, filenames)
        write_multiple_groups_reads(A, B, filenames)
        allReadsID = write_ambiguous_reads(C, filenames)
        D = allReadsID - C
        for i in A:
            D = D - A[i]
            D = D - B[i]

        write_nonSNV_reads(D, filenames)
    elif filenames[0].endswith("fasta"):
        #####write_two_groups_fasta(A, B, filenames)
        write_multiple_groups_fasta(A, B, filenames)
        allReadsID =write_ambiguous_fasta(C, filenames)
        D = allReadsID - C
        for i in A:
            D = D - A[i]
            D = D - B[i]

        write_nonSNV_fasta(D, filenames)
    '''
        

#-------useless-------

def try_merge_block_according_TGS_reads(readID, nodeID, groups1, groups2):

    size = len(groups1)

    inter = []
    for i in range(size):
        l1 = set( nodeID ).intersection( set(groups1[i]) ) 
        l2 = set( nodeID ).intersection( set(groups2[i]) )
        if len(l1) > 0:
            inter.append(str(i+1)+'_A_'+str(len(l1)))
        if len(l2) > 0:
            inter.append(str(i+1)+'_B_'+str(len(l2)))
    if len(inter) > 2 or (len(inter)==2 and inter[0].split('_')[0] != inter[1].split('_')[0] ):
        print (readID, inter, nodeID)
    return

