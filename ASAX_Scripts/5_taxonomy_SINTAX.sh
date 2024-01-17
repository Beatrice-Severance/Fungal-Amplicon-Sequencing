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
vsearch -sintax /scratch/aubbxs/GAUSDA/2023/output/cluster/otus_R1.fasta -db /home/aubbxs/noel_work/EVSMITH/fungi/tax_databases/utax_reference_dataset_SINTAX_29.11.2022_all_Euk_SynMock.fasta -tabbedout /scratch/aubbxs/GAUSDA/2023/output/SINTAX/fungi_R1_UNITE_new.txt -strand both -sintax_cutoff 0.8
