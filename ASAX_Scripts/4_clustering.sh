#!/bin/bash

############################################
#     Bioinformatic pipeline     	   #
#     ITS amplicon sequences	           #
#DEREPLICATION, CLUSTERING, CHIMERA REMOVAL#
############################################

module load vsearch

# dereplication 
vsearch --derep_fulllength /scratch/aubbxs/GAUSDA/2023/output/filtertrim/trimmed_R1.fasta --output /scratch/aubbxs/GAUSDA/2023/output/cluster/uniques_R1.fasta -sizeout

# de-noising (error correction), output is zero-radius OTUs
usearch -unoise3 /scratch/aubbxs/GAUSDA/2023/output/cluster/uniques_R1.fasta -tabbedout /scratch/aubbxs/GAUSDA/2023/output/cluster/unoise_zotus_R1.txt -zotus /scratch/aubbxs/GAUSDA/2023/output/cluster/zotus_R1.fasta

# clusters OTUs based on traditional 97% identity 
usearch -cluster_otus /scratch/aubbxs/GAUSDA/2023/output/cluster/uniques_R1.fasta -minsize 2 -otus /scratch/aubbxs/GAUSDA/2023/output/cluster/otus_R1.fasta -uparseout /scratch/aubbxs/GAUSDA/2023/output/cluster/uparse_otus_R1.txt -relabel FOTU_ --threads 20
