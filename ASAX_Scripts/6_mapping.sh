#!/bin/bash 

############################################
#     Bioinformatic pipeline               #
#     ITS amplicon sequences               #
#            MAPPING                       #
############################################

# First use seqtk to convert all the demultiplexed samples into fasta on a loop.
# Use samples.txt like we did for cutadapt

# Load the modules
module load gcc/11.2.0
module load libiconv/1.16
module load plumed/2.6.3
module load velvet/1.2.10
module load seqtk/1.3-olt7cls
source /apps/profiles/modules_asax.sh.dyn
module load python/3.3.2
module load anaconda/3-2023.03

for sample in $(cat samples.txt)
do

echo "On sample: $sample"
    seqtk seq -a /scratch/aubbxs/GAUSDA/2023/fwreads/${sample}*.fastq.gz > ${sample}.fasta

    # have to replace the beginning of the fasta headers with the file name for mapping. Otherwise we # get one sample with all the read counts, which is not what we want.

    python3 replacefastaheaders_filename.py ${sample}.fasta

done

# have to create one file containing all the reads from the demultiplexed reads
cat *_newheaders.fasta > demultiplexed_new.fasta


# align the demultiplexed reads back to the now clustered OTUs or ZOTUs (ESV)
module load anaconda/3-2021.11
module load vsearch
module load usearch/11.0.667

vsearch -usearch_global demultiplexed_new.fasta -db /scratch/aubbxs/GAUSDA/2023/output/cluster/otus_R1.fasta -strand plus -id 0.97 -otutabout otu_table_ITS_UNOISE_R1.txt
