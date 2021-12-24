########################################################################
# File Name: ../step4_mutiple_thread.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Tue 14 May 2019 16:36:45 AEST
#########################################################################
#!/bin/bash
import threading
import sys
##import Queue
import queue
import tools
from multiprocessing import Process
# https://blog.csdn.net/zx8167107/article/details/81083249

#import build_matrix_for_PairEndReads


def get_longest_consist_len( l ):
    # [1, 2, 3, 4, 5, 6, 7, 8, 20, 21, 22]
    # return 8, 1
    sortedL = sorted(l)
    lenL = len(l)
    cnt, maxCnt = 1, 1
    pre, start, maxStart = l[0], l[0], l[0]
    #print (l)
    for i in range(1,lenL):
        if l[i] == pre+1:
            cnt += 1
        else:
            start = l[i]
            cnt = 1
        pre = l[i] 
        if cnt > maxCnt:
            maxCnt = cnt
            maxStart = start
    #print (maxCnt, maxStart)    
    return maxCnt, maxStart





#class Reader(threading.Thread):
class Reader(Process):
    def __init__(self, thread_id, NGS, temp_list, temp_list2, readType=1, k=21):
        super(Reader, self).__init__()
        self.thread_id = thread_id
        self.NGS = NGS
        self.temp_list = temp_list
        #self.temp_list2 = temp_list2
        self.K = k
        self.readType = readType

    def scan_reads_search_kmer(self, seq):

        seqLen = len(seq)
        intersection = {} # []
        pre1, pre2 = (-1, -1), (-1, -1)
        for i in range(seqLen - self.K):
            key = str(seq[i:i + self.K])
            if key.count('A') + key.count('C') + key.count('T') + key.count('G')  != self.K:
                print (key, "reads include other elements (except ATCG)")
                continue
            Rkey = tools.reverse(key)
            Rflag = False 
            # reduce one binarySearch
            if key > Rkey:
                key = Rkey
                Rflag = True

            re1 = tools.binarySearch(self.NGS, key)

            ''' # version one
            if re1 != -1: 
                #intersection.append( (re1[1], re1[2]+"f") )
                if re1[1] != pre1[1]:
                    if Rflag == False:
                        intersection.append( (re1[1], "f") ) # f and r here no meaning now
                    else:
                        intersection.append( (re1[1], "r") )
                    pre1 = re1
            '''
            # version two, at least two happen
            if re1 != -1: 
                if re1[1] not in intersection:
                    intersection[ re1[1] ] = []
                intersection[ re1[1] ].append( i ) 

        inter = []   
        if len(intersection) > 0:
            for key in intersection:
                #print (key, intersection[key])
                #supportPosNum =  len(intersection[key]) #if supportPosNum > 31: #maybe not consis position
                supportPosNum, startPos = get_longest_consist_len( intersection[key] ) #if supportPosNum > 31: #maybe not consis position
                inter.append( (key, str(supportPosNum) + '_' + str(startPos) ) )
                #if supportPosNum == 1:
                    #sys.exit()

                # heterozygous kmer id, and heterozygous binay code, pos on reads
            #re2 = binarySearch(self.NGS, Rkey)
            #if re2 != -1: #print "pos:", i, re2
                #intersection.append( (re2[1], re2[2]+"r") )
                #if re2[1] != pre2[1]:
                    #intersection.append( (re2[1], "r") )
                    #pre2 = re2
        #else:
            #print ("no intersection, seq length", len(seq))
        return inter


    def run(self):
        state = 0
        fout = open("part_matrix_"+ str(self.thread_id), "w")
        count = 0
        textLen = len(self.temp_list)
        for i in range(textLen):
            text = self.temp_list[i]
            if state==0 and (text.startswith(">") or text.startswith("@") ):
                readID = text.strip()[1:].split(' ')[0]
                state = 1
                count += 1
                if count % 1000 == 0:
                    print ("thread ", self.thread_id, "deal reads ", count) 
            elif state == 1:
                state = 0
                seq = text.strip()
                #print (readID)
                intersection = self.scan_reads_search_kmer(seq)
                if len(intersection) < 1:
                    print ("%s no hetezygous kmers on this read, read length %s" % (readID, len(seq) ) )
                    continue
                fout.write( "%s %s " % ( len(intersection), readID ) )
                for (p, binary) in intersection:
                    fout.write( "%s %s " % (p, binary) )
                score=len(intersection)*'4'    
                fout.write("%s\n" % score)
            else:
                #print (text)
                # for fasta, no else; for fastq seq score
                continue
        fout.close()
 
 
class Partition(object):
    def __init__(self, file_name, thread_num):
        self.file_name = file_name
        self.block_num = thread_num
 
    def part_and_queue(self):
        pos_list = []
        file_size = tools.file_lines(self.file_name)
        block_size = int(file_size / self.block_num) # each reads seq should in one line
        if block_size % 2 == 1: # make sure read.seq and readID in same bin 
            block_size += 1
        start_pos = 0
        #global q 
        for i in range(self.block_num):
            if i == self.block_num - 1:
                end_pos = file_size - 1
                pos_list.append((start_pos, end_pos))
                break
            end_pos = start_pos + block_size - 1
            if end_pos >= file_size:
                end_pos = file_size - 1
            if start_pos >= file_size:
                break
            pos_list.append((start_pos, end_pos))
            start_pos = end_pos + 1
 
        print ("partition list", pos_list)
        fd = open(self.file_name, 'r')
        for pos_tu in pos_list:
            temp_text = []
            start = pos_tu[0]
            end = pos_tu[1]
            while start <= end:
                text = fd.readline().strip('\n')
                temp_text.append(text)
                start = start + 1
            q.put( (temp_text, "") )
        print ("q size", q.qsize())    
        fd.close()
 


def read_kmerID(filename): # for extend kmer, list

    heteKmer = []
    with open(filename, "r") as f:
        for line in f:
            heteKmer.append( line.strip().split() )
    return sorted(heteKmer) # should not have same kmer, diff ID
       

def build(kmerFile, TGS_files, thread_num, readType, k):

    TGS_filenames = TGS_files.split(',')
    TGS_filesNum = len(TGS_filenames)
    heteKmer = read_kmerID(kmerFile)
    if TGS_filesNum == 0:
        print ("long reads input error")
        sys.exit()

    print ( "heterozgyous kmer number:",len(heteKmer) ) 
    #thread number
    global q
    q = queue.Queue()
    #p = Partition(pair_filename1, pair_filename2, thread_num)
    for f in TGS_filenames:
        #assert thread_num >= TGS_filesNum
        p = Partition(f, int(thread_num/TGS_filesNum))
        p.part_and_queue()
    
    t = []
    for i in range(thread_num):
        (temp_list, temp_list2) = q.get()
        t.append( Reader(i, heteKmer, temp_list, temp_list2, readType, k) )
    for i in range(thread_num):
        t[i].start()
    for i in range(thread_num):
        t[i].join()


if __name__ == '__main__':
    
    if len(sys.argv) < 5:
        print ("chr22.NGS.kmer long_insertion_NGS(TGS).fastq thred_num readType(default 1) k")
        sys.exit()

    NGS_filename = sys.argv[1]
    TGS_filenames = sys.argv[2].split(',')
    TGS_filesNum = len(TGS_filenames)
    readType = int(sys.argv[4])
    k = int(sys.argv[5])
    
    #filesNum = sys.argv[2].count(',') + 1
    #for i in range(filesNum):
    #    TGS_filename = sys.argv[2]
    #pair_filename2 = ""
    #heteKmer = build_matrix_for_PairEndReads.kmerPair2KmerID(NGS_filename)
    heteKmer = read_kmerID(NGS_filename)
    print ( "heterozgyous kmer number:",len(heteKmer) ) 
 
    #thread number
    thread_num = int(sys.argv[3]) 
    global q
    q = queue.Queue()
    #p = Partition(pair_filename1, pair_filename2, thread_num)
    for f in TGS_filenames:
        #assert thread_num >= TGS_filesNum
        p = Partition(f, int(thread_num/TGS_filesNum))
        p.part_and_queue()
    
    t = []
    for i in range(thread_num):
        (temp_list, temp_list2) = q.get()
        t.append( Reader(i, heteKmer, temp_list, temp_list2, readType, k) )
    for i in range(thread_num):
        t[i].start()
    for i in range(thread_num):
        t[i].join()


