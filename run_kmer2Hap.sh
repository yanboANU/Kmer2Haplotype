#########################################################################
# File Name: runkmercalling.sh
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Thu 04 Jul 2019 16:52:18 AEST
#########################################################################
#!/bin/bash
set -e
#cd $$3x/chr$1/
#mkdir kmercalling
#cd kmercalling

if [ $# != 8 ]; then
	echo "\$1:chrID \$2:k \$3:k-1 \$4:hom cov \$5:NGS.fq(fa) \$6:TGS.fq \$7 genomeSize \$8:maxThreads" 
	exit 
fi	

echo "\$1:chrID \$2:k \$3:k-1 \$4:hom cov \$5:NGS.fq(fa) \$6:TGS.fq \$7 genomeSize \$8:maxThreads" 
echo $1 $2 $3 $4 $5 $6 $7 $8 

species=$1
TGS=$6
genomeSize=$7
maxThreads=$8

start=`date +%s`

#if [ ! -f "k_31_pair.snp" ]; then
#	sh /home/yanbo/software/myVersionKmer2SNP/runkmercalling.sh $species $2 $3 $4 $5 0
#else
#        echo "kmer2snp already done"
#fi
foldername="./Kmer2Hap"

if [ ! -d $foldername  ]; then
	mkdir $foldername
else
	echo foldername, "folder exist"
fi	

cd $foldername
pwd

path="/home/yanbo/software/Kmer2Haplotype"
if [ ! -f "all_matrix" ]; then
    /home/yulin/py36/bin/python3 $path/libprism/local/using_extend_kmer.py ../../k_31_pair.snp ../../k_31_pair.snp.extend

    /home/yulin/py36/bin/python3 $path/libprism/local/build_matrix_for_supperReads.py kmer_ID $TGS $maxThreads 1 31 >step1.log

    cat par* >all_matrix
	#rm par*_matrix
else
	echo "already genenrate matrix"
fi

if [ ! -f "phased_kmer_ID" ]; then
    /home/yulin/py36/bin/python3 $path/libprism/local/build_graph_notrick.py all_matrix 100 >build_graph.log
fi 

if [ ! -f "groupA.fastq" ]; then
    /home/yulin/py36/bin/python3 $path/libprism/local/binning_TGS_reads.py phased_kmer_ID all_matrix $TGS > binning.log
fi 

foldername="./canu_assembly"

if [ ! -d $foldername  ]; then
	mkdir $foldername
else
	echo foldername, "folder exist"
fi	
cd $foldername

pwd

sh $path/script/run_canu.sh $species $genomeSize $maxThreads 

end=`date +%s`
runtime=$((end-start))
echo run $runtime seconds



#/home/yulin/py36/bin/python3 /home/yanbo/software/Kmer2Haplotype/convert_output.py kmer_ID_mid phased_kmer_ID >convert_output.log

#python3 /home/yanbo/software/Kmer2Haplotype/libprism/evaluate/compare.py fungi phased_kmer >compare.log
