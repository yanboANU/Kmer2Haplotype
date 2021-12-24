#########################################################################
# File Name: tools.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 09 May 2019 13:12:45 AEST
#########################################################################
#!/bin/bash

import sys
import os
#import subprocess
import networkx as nx

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import threading
import queue

class MyThread(threading.Thread):
    def __init__(self, func, args=()):
        #super(MyThread,self).__int__()
        threading.Thread.__init__(self)
        self.func = func
        self.args = args
    def run(self):
        #self.result = self.func(*self.args)
        self.result = self.func(self.args)
    def get_result(self):
        try:
            return self.result
        except Exception:
            return None
 

m={}
m['A'] = 'T'
m['T'] = 'A'
m['C'] = 'G'
m['G'] = 'C'

def reverse(s): # AATC
    news = ""
    for ele in s:
        news = m[ele] + news #CGATT
    return news   


def reverse2(s): # AATCG
    news = ""
    #print ("total length",len(s))
    #count = 0
    for ele in s:
       if ele in m:
           news = m[ele] + news #CGATT
       else:
           news = ele + news #N,M, in fungi reads  
       #count += 1
       #if count % 100000 == 0:
           #print (count)
    return news  


def reverse3(s): # reverse2 and reverse3 cost similar time
    news = ""
    for ele in iter(s):
       if ele ==  'A':
           news = 'T' + news #CGATT
       elif ele == 'T':
           news = 'A' + news 
       elif ele == 'C':
           news = 'G' + news 
       elif ele == 'G':
           news = 'C' + news 
       else: #N,M, in fungi reads 
           news = ele + news
    return news  


def reverse_long_string(s): # mush faster than reverse2 and reverse3
    l = len(s)
    if l <= 10000:
        return reverse3(s)

    thread_num = min(30, int(l/10000)) 
    global q
    q = queue.Queue()
    p = Partition(s, l, thread_num)
    p.part_and_queue()
    
    tList = []
    for i in range(thread_num):
        (s) = q.get()

        #print ("check", type(s))
        #print ("check", s)
        #sys.exit()
        t = MyThread(reverse3, args=(str(s)))
        tList.append( t )

    for i in range(thread_num):
        tList[i].start()
     
    res = {} 
    for i in range(thread_num):
        tList[i].join()
        res[i] = tList[i].get_result()
    final = ""
    for i in range(thread_num):
        final = res[i] + final
    return final    


class Partition(object):
    def __init__(self, s, l, thread_num):
        self.s = s
        self.l = l
        self.thread_num = thread_num
 
    def part_and_queue(self):
        pos_list = []
        block_size = int(self.l / self.thread_num)  # each reads seq should in one line
        start_pos = 0
        #global q
 
        for i in range(self.thread_num):
            if i == self.thread_num - 1:
                end_pos = self.l #- 1
                pos_list.append((start_pos, end_pos))
                break
            end_pos = start_pos + block_size - 1
            if end_pos >= self.l:
                end_pos = self.l #- 1
            if start_pos >= self.l:
                break
            pos_list.append((start_pos, end_pos))
            start_pos = end_pos # + 1
 
        for pos_tu in pos_list:
            start = pos_tu[0]
            end = pos_tu[1]
            temp_text = self.s[start:end]
            #print ("check", type(temp_text))
            #print ("check", temp_text)
            #sys.exit()
            q.put( temp_text )
        #print ("q size", q.qsize())    
 




def reverse_ward(ward):
    if ward == 'f':
        return 'b'
    elif ward == 'b':
        return 'f'
    else:
        print ("reverse forward or backward error")
        sys.exit()



def get_smaller_pair_kmer_keep_order(k1, k2):
    

    Rk1 = reverse(k1)
    Rk2 = reverse(k2)
    if (Rk1 < k1 and Rk1 < k2) : # the mini is Rk1 or Rk2
        k1,k2 = Rk1, Rk2    
    if (Rk2 < k1 and Rk2 < k2) : # the mini is Rk1 or Rk2
        k1,k2 = Rk1, Rk2
    return k1, k2

def get_smaller_pair_kmer(kmer1, kmer2):
   
    # kmer1, and kmer2 may switch
    lenKmer1 = len(kmer1)
    lenKmer2 = len(kmer2) 
    if lenKmer1 < lenKmer2:
        kmer1, kmer2 = kmer2, kmer1 
    RKmer1 = reverse(kmer1)
    RKmer2 = reverse(kmer2)
    if lenKmer1 == lenKmer2:
        if kmer2 < kmer1:
            kmer1, kmer2 = kmer2, kmer1
        if RKmer2 < RKmer1:
            RKmer1, RKmer2 = RKmer2, RKmer1
    if RKmer1 < kmer1:
        kmer1, kmer2 = RKmer1, RKmer2
    return kmer1, kmer2               
# above code from myversionKmer2SNP

def get_smaller_kmer(kmer):
   
    RKmer = reverse(kmer)
    if kmer < RKmer:
        return kmer
    else:
        return RKmer


################################
#fasta and fastq
###############################
def read_fasta_file(filename):

    seqs = {}
    for seq_re in SeqIO.parse(filename, "fasta"):
        seqs[ seq_re.id ] = seq_re.seq
        #print (seq_re.id)
    return seqs

def print_seq(seqID, startPos, length, seqs):
    print ("seqID: %s, seq length: %s, start position: %s" % (seqID, len(seqs[seqID]), startPos) )
    #fseq = seqs[seqID][startPos : startPos + length].upper()
    fseq = seqs[seqID][startPos - length : startPos + length].upper()
    print ( fseq )
    print ( reverse( fseq ) )

    #bseq = ( reverse( seqs[seqID].upper() ) )[startPos : startPos + length].upper()
    #print (bseq)
    #print ( reverse( bseq ) )
    
    return 



def write_one_fastq_record(fout, seqRecord):
    
    fout.write("@%s\n" % seqRecord.id)
    fout.write("%s\n" %seqRecord.seq)
    fout.write("+%s\n" % seqRecord.id)
    fout.write("%s\n" % seqRecord.format("fastq").split('\n')[-2] )
    return


def write_one_fasta_record(fout, seqRecord):
    
    fout.write(">%s\n" % seqRecord.id)
    fout.write("%s\n" %seqRecord.seq)
    return


def write_list(fout, l):
    #print (l)
    for ele in l[:-1]:
        fout.write("%s " % ele)
    fout.write("%s\n" % l[len(l)-1])
    return

def check_node(G, node):
   
    node1= node + "_0"
    print (node1, nx.degree(G, node1 )  )
    for ele in list(G.neighbors(node1) ):
            w = G.get_edge_data(node1, ele)['weight']
            #print (ele, w, end=" ")
    print ("")
    node2= node + "_1"
    print (node2, nx.degree(G, node2)  )
    for ele in list(G.neighbors(node2) ):
            w = G.get_edge_data(node2, ele)['weight']
            #print (ele, w, end=" ")
    print ("")  

def update_dict(cnt, n):
    
    if n not in cnt:
        cnt[n] = 1
    else:
        cnt[n] += 1


def check_groups_mix(group, n1):

    chrID = set()
    if isinstance(group[n1], set):
        for node in group[n1]:
            if node.count('_')==2:
                chrID.add(node[:2])          
    elif isinstance(group[ group[n1] ], set):
        for node in group[ group[n1] ]:
            if node.count('_')==2:
                chrID.add(node[:2])
    if len(chrID) > 1:
        return True
    else:
        return False


def print_final_result(finalGroups):
    
    print ("final group:", len(finalGroups) )
    cnt = 0 
    for g in finalGroups:
        cnt += len(g)
        print (g)
    print ("node number in final group:", cnt )
    return

def assign_ID_2_kmerPair(kmerPair):
    
    heteKmer = []
    count = 1 # hete kmer ID
    for (k1, k2, w) in kmerPair:
        heteKmer.append( (k1, count, '0') ) # (k1, str(count) + '_0')
        heteKmer.append( (k2, count, '1') )
        count += 1
    return sorted(heteKmer)




def binarySearch(array, t):
    low = 0
    height = len(array)-1
    while low <= height: # debug should be <=
        mid = low + int( (height-low)/2 )
        if array[mid][0] < t:
            low = mid + 1
        elif array[mid][0] > t:
            height = mid - 1
        else:
            return array[mid]
    return -1

def binarySearch_ListValue(array, t):
    # return first index, value > t
    low = 0
    height = len(array)-1
    while low <= height: 
        mid = low + int( (height-low)/2 )
        if len(array[mid][1]) <= t:
            low = mid + 1
        elif len(array[mid][1]) > t:
            height = mid - 1
    if low >= len(array):
        return low
    elif len(array[low][1])==t:
        low += 1
    return low




'''
s="ATCG"
print s
print reverse(s)
'''
 

def check_node_set(nodeSet):

    newSet = set()
    for n in nodeSet:
        symN = get_symmetrical_node(n)
        if symN not in nodeSet:
            newSet.add(symN)

    return sorted( newSet.union(nodeSet) )
    

def is_symmetrical_node(n1, n2):
    if n1[:-1] == n2[:-1] and n1[-1] != n2[-1]:
        return True
    return False

def get_symmetrical_node(node):
    if node.endswith('1'):
        return node[:-1] + '0'

    if node.endswith('0'):
        return node[:-1] + '1'

def print_list(l):
    for ele in l:
        print(ele,)
    print ("")

def file_lines(filename):
    command = "wc -l " + filename + " >file_lines"
    os.system(command)
    with open("file_lines", "r") as f:
        for line in f:
            words = line.split()
            return int(words[0])


def hamming_distance(s1, s2):
    count = 0
    lenS = len(s1)
    #print s1, s2
    assert lenS == len(s2)
    for i in range(lenS):
        if s1[i] == -1 or s2[i] == -1:
            continue
        if s1[i] != s2[i]:
            count +=1
    return count

def merge_pair_end_reads(r1, r2, pairEndType):
    
    # pair end reads with overlap 
    #        <- r1
    # r2 ->
    # reads length 150
    # insert size 200, variation 10
    seq = r2
    if pairEndType == 0: # pair end reads with overlap
        overlapLen = 100
        Rr1 = reverse(r1)  
        overlap = r2[-overlapLen :] 
        pos = Rr1.find(overlap)
        while pos == -1:
            overlapLen -= 1
            overlap = r2[-overlapLen :] 
            pos = Rr1.find(overlap)
        if pos > 0:
            seq = r2 + Rr1[overlapLen + pos : ]
            #print r2 #print Rr1 #print seq
    if pairEndType == 1:  # mate pair reads, no overlap
        seq = reverse(r1) + r2 # 
    return seq    

def get_non_zero_index(l): # list l
    
    index = []
    for i in range(len(l)):
        if l[i] != 0:
            index.append(i)
    return index        
