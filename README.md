# Kmer2Haplotype
**Kmer2Haplotype** is a pipeline of reference-free haplotype phasing by incorporating both short and long reads from the same individual.

# Downloading **Kmer2Haplotype**
To download Kmer2SNP, you have to install the following software:


<pre><code>
(1) Kmer2SNP: 
https://github.com/yanboANU/Kmer2SNP
(2) Flye:
https://github.com/fenderglass/Flye
(3) python3
 </code></pre>
  
Then clone the **Kmer2Haplotype** repository to your machine.
<pre><code> git clone https://github.com/yanboANU/Kmer2Haplotype.git </code></pre>  

# Input Format
The input of **Kmer2Haplotype**: 
<pre><code>
(1) NGS short reads: can be fasta, fastq, either gzipped or not.  
(2) TGS long reads: can be fasta, fastq
</code></pre>


# Output Format
<pre><code>
H1_contigs.fasta
H2_contigs.fasta
unphasd.fasta
</code></pre>

# Example Usage
<pre><code>
python3 kmer2Hap.py
</code></pre>
