#!/bin/bash -login

# STEP 5 - MAPPING

############################################
#     Bioinformatic pipeline               #
#     ITS amplicon sequences               #
#            MAPPING                       #
############################################

# First use seqtk to convert all the demultiplexed samples into fasta on a loop.
# Use samples.txt like we did for cutadapt

# Load the modules
module load gcc/11.2.0
module load seqtk/1.3-olt7cls

for sample in $(cat combsamples.txt)
do

echo "On sample: $sample"
    seqtk seq -a /scratch/aubbxs/fungi/21-22/21-22forward/${sample}*.fastq.gz > ${sample}.fasta

    # have to replace the beginning of the fasta headers with the file name for mapping. Otherwise we get one sample with all the read counts, which is not what we want.
    # We use awk to append the filename at the beginning of each fasta sequence after the >, then we pipe it to sed to replace the underscore with a period.

    awk '/>/{sub(">","&"FILENAME":");sub(/\.fasta/,x)}1' ${sample}.fasta | sed '/^>/s/_/\ /g' > ${sample}_new.fasta

done

# have to create one file containing all the reads from the demultiplexed reads
cat *_new.fasta > demultiplexed_new.fasta

#move demultiplexed_new.fasta to mapped folder
mv *demultiplexed_new.fasta /scratch/aubbxs/fungi/21-22/output/map/

# Taking out the /output/demultiplexed/ text from the beginning of the sample names - probably a better way to run this, but it works. 
#sed 's/[/home/aublxr001/p.soil_fungi/Demultiplexed/]//g' /home/aublxr001/p.soil_fungi/Demultiplexed/demultiplexed_new.fasta > /home/aublxr001/p.soil_fungi/Demultiplexed/#demultiplexed_new2.fasta

awk '{gsub("/scratch/aubbxs/fungi/21-22/output/map/",""); print}' /scratch/aubbxs/fungi/21-22/output/map/demultiplexed_new.fasta | sed '/^>/s/_/\ /g' > demultiplexed_new2.fasta

# align the demultiplexed reads back to the now clustered OTUs or ZOTUs (ESV)
module load vsearch
vsearch -usearch_global demultiplexed_new2.fasta -db /scratch/aubbxs/fungi/21-22/output/cluster/otus_R1.fasta -strand plus -id 0.97 -otutabout otu_table_ITS_UNOISE_R1.txt
