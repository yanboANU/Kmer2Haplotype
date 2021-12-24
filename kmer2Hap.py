#########################################################################
# File Name: kmer2Hap.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Wed 15 Jan 2020 09:56:19 AEDT
#########################################################################
#!/bin/bash
import os
import sys


import argparse
import os
import sys
import logging
import time
from libprism.local import process
#import networkx as nx
#import draw
#############################################################################


def run(args):

    command = "sh /g/data/xl04/yl4079/code/Kmer2Hap/script/run_kmer2snp.sh " + args.k + " " + args.c1 + " " + args.NGS
    print (command)
    os.system( command )
    if os.path.exists("Kmer2Hap")== False:
        os.system( "mkdir Kmer2Hap")
    os.chdir("Kmer2Hap")

    if args.n == None:
        args.n = 10
    if args.m == None:
        args.m = 1

    if args.real == None:
        args.real = 1
    #print ("debug real", args.real)
    process.run(args.k, args.TGS, args.n, args.m, args.c2, bool(int(args.real)))

    pwd = os.getcwd()
    path = pwd + "/group_fastq"
    os.chdir(path)
    
    #sys.exit()
    '''
    if os.path.exists("canu_assembly")== False:
        os.system( "mkdir canu_assembly")
    os.chdir("canu_assembly")

    #command = "sh /home/yanbo/software/Kmer2Haplotype/script/run_canu.sh "  + path
    command  = "sh /g/data/xl04/yl4079/code/Kmer2Hap/script/run_canu.sh "  + path
    '''

    if os.path.exists("flye_assembly")== False:
        os.system( "mkdir flye_assembly")
    os.chdir("flye_assembly")
    #command = "sh /home/yanbo/software/Kmer2Haplotype/script/run_flye.sh "  + path
    command  = "sh /g/data/xl04/yl4079/code/Kmer2Hap/script/run_flye.sh "  + path
    # + args.r + " " +  args.l + " "  + str(args.n)
    print (command)
    os.system( command )


parser = argparse.ArgumentParser(description="Haplotype assembly for Single-individual based on NGS reads")

parser.add_argument('--r', help='chromosome ID or sampleID', required=True)
parser.add_argument('--k', help='kmer size', required=True)
parser.add_argument('--NGS', help='next-generation sequencing (NGS) reads', required=True)
parser.add_argument('--c1', help='homozygous coverage of NGS reads', required=True)
parser.add_argument('--TGS', help='third generation sequencing(TGS) reads', required=True)
parser.add_argument('--c2', help='homozygous coverage of TGS reads', required=True)
parser.add_argument('--l', help='genome Size', required=True)
parser.add_argument('--real', help='when realData set with lower NGS read coverage, set True (Default True)')
parser.add_argument('--n', help='maximum threads (default 10)', required=False)
parser.add_argument('--m', help='matrix type, for m=0, keep all infor, m=1, remove some low coverage reads (default 1)', required=False)
parser.add_argument('--output_dir', help='output directory (Default is ./)', required=False)
#parser.add_argument('--o', help='output file prefix', required=True)
#parser.add_argument('--ref', help='reference or scaffolds', required=True)



args = parser.parse_args()
print ('Kmer2Hap Input Arguments:', ','.join([ '{}={}'.format(i, vars(args)[i]) for i in vars(args)]))


#if __name__ == '__main__':
    
if args.output_dir == None:
    args.output_dir = './'
os.chdir( args.output_dir )
run(args)
    



