#1_strippingPrimers.sh
#!/bin/bash
#


##################################
#     Bioinformatic pipeline     #
#     ITS amplicon sequences     #
#                                #
#     ----------------------     #
#        stripping primers       #
#                                #
##################################

# This script is run from the same directory as the sequence files in the scratch. It loops over the names in samples.txt
# samples.txt is made by simply running the command ls > samples.txt then using nano and editing the last listed file since this is samples.txt

#  load the module
module load anaconda/2-4.2.0_cent

for sample in $(cat 2022samples.txt)
do

    echo "On sample: $sample"

    cutadapt -g CTTGGTCATTTAGAGGAAGTAA -e 0.01 --discard-untrimmed --match-read-wildcards /scratch/aubbxs/fungi/Demultiplexed_QUO1003506_PecanNCST_Fungi/forward_reads/${sample}*.fastq.gz >> ${sample}_trimmed.fastq

done

# -n trimming more than one adapter each read
# -e expected error rate (default 10%)
# --discard-untrimmed discards the pair if one of the reads does not contain an adapter
# --match-read-wildcards All IUPAC nucleotide codes (wildcard characters) are supported
