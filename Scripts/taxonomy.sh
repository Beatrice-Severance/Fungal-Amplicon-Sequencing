#5_DADA2_NBC.sh
####RUN THIS SCRIPT ON LARGE WITH 1 CORE AND 32GB OF MEMORY! It is a massive memory hog and the job will kill itself if you don't give around this much
#!/bin/bash

############################################
#     Bioinformatic pipeline               #
#     R analysis ANCOM                     #
#   Differential abundance test            #
# To run this script you must edit the ancom_v2.1.R script to the appropriate variables    #
############################################

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load R/4.1.0

R CMD BATCH dada2_assigntax_NBC.R
