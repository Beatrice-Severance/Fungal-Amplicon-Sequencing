#!/bin/bash

##################################
#     Bioinformatic pipeline     #
#     ITS FUNGI 		 #
#     ----------------------     #
#     taxonomy assignment	 #
#                                #
#      Zachary Noel      	 #
#    Michigan State University   #
##################################

module load vsearch

# Assign taxonomy using SINTAX algorithm
vsearch -sintax /scratch/aubbxs/GAUSDA/2023/output/cluster/otus_R1.fasta -db /home/aubbxs/noel_work/EVSMITH/fungi/tax_databases/utax_reference_dataset_SINTAX_29.11.2022_all_Euk_SynMock.fasta -tabbedout fungi_R1_UNITE.txt -strand both -sintax_cutoff 0.8

# Print the first and fourth columns of the sintax output, replace the commas with tabs, get rid of the taxonomic prefixes, save it in a text file
awk '{print $1, $4}' fungi_R1_UNITE.txt | tr ' ' ',' |  sed 's/d://; s/p://; s/c://; s/o://; s/f://; s/g://; s/s://' > fungi_R1_UNITE2.csv
