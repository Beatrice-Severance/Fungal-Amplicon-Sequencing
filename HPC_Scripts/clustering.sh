#4_clustering.sh
#!/bin/bash

############################################
#     Bioinformatic pipeline               #
#     ITS amplicon sequences               #
#DEREPLICATION, CLUSTERING, CHIMERA REMOVAL#
############################################

#  HAVE TO INSTALL USEARCH BEFORE USING    #

module load vsearch

# dereplication
vsearch --derep_fulllength /scratch/aubbxs/fungi/output/filtered/trimmed_R1.fasta --output /scratch/aubbxs/fungi/output/clustered/uniques_R1.fasta -sizeout

# de-noising (error correction), output is zero-radius OTUs
usearch -unoise3 /scratch/aubbxs/fungi/output/clustered/uniques_R1.fasta -tabbedout /scratch/aubbxs/fungi/output/clustered/unoise_zotus_R1.txt -zotus /scratch/aubbxs/fungi/output/clustered/zotus_R1.fasta

# clusters OTUs based on traditional 97% identity
usearch -cluster_otus /scratch/aubbxs/fungi/output/clustered/uniques_R1.fasta -minsize 2 -otus /scratch/aubbxs/fungi/output/clustered/otus_R1.fasta -uparseout /scratch/aubbxs/fungi/output/clustered/uparse_otus_R1.txt -relabel FOTU_ --threads 20
