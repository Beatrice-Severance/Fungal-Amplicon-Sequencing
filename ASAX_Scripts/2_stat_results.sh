#!/bin/bash

############################################
#     Bioinformatic pipeline     	   #
#     ITS amplicon sequences	           #
#	  	STATS  		  	   #
############################################

# Usually just trim at 250 bp and maxee 1.0

cat /scratch/aubbxs/GAUSDA/output/strip/*trimmed.fastq > /scratch/aubbxs/GAUSDA/output/stats/trimmed_combined.fastq

module load anaconda/3-2021.11

vsearch -fastq_stats /scratch/aubbxs/GAUSDA/output/stats/trimmed_combined.fastq -log /scratch/aubbxs/GAUSDA/output/stats/stats_results_R1.txt
