#########################################################################
# File Name: libprism/local/process.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Mon 28 Jun 2021 06:10:18 PM AEST
#########################################################################
#!/bin/bash
import os
import sys

from libprism.local import using_extend_kmer
from libprism.local import build_matrix_for_supperReads
from libprism.local import build_graph_version2
from libprism.local import binning_TGS_reads




def run(k, longReadFiles, threadNum, mType, c2, l):

    print ("in process run, input datatype(1 real; 0 simulated) ", l) 
    k = int(k)
    threadNum = int(threadNum)
    if os.path.exists("kmer_ID") == False:
        file1 = "../k_" + str(k) + "_pair.snp"
        file2 = "../k_" + str(k) + "_pair.snp.extend"
        using_extend_kmer.get_heterozygous_snp_kmer(file1, file2)
    
    if os.path.exists("all_matrix") == False:
        file1 = "kmer_ID"
        readType = 1 # 1 means TGS reads, now only support TGS reads 
        build_matrix_for_supperReads.build(file1, longReadFiles , threadNum, readType, k)
        os.system("cat par* >all_matrix")
        os.system("rm par*_matrix*")
    
    if os.path.exists("phased_kmer_ID") == False:    
        #/home/yulin/py36/bin/python3 $path/libprism/local/build_graph_version2.py all_matrix $k $filterRate $homCov 0 >build_graph.log
        mType = int(mType)
        c2 = int(c2)

        #print ("debug in run", l) 
        build_graph_version2.phasing_kmer(mType, k, c2, l)
  
    #sys.exit() 
    if os.path.exists("group_fastq") == False and os.path.exists("group_fasta") == False :
        print ("binning TGS reads")
        binning_TGS_reads.run("phased_kmer_ID", "all_matrix", longReadFiles, k)
    
