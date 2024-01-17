#2_stat_results.sh
#!/bin/bash

############################################
#     Bioinformatic pipeline               #
#     ITS amplicon sequences               #
#               STATS                      #
############################################

# Usually just trim at 250 bp and maxee 1.0

cat /scratch/aubbxs/fungi/output/trimmed/*trimmed.fastq > trimmed_combined.fastq

module load vsearch
vsearch -fastq_stats trimmed_combined.fastq -log stats_results_R1.txt
