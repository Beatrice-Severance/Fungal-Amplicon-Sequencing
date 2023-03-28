#run the below two lines before running the script
#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load R/4.1.0

library(dada2)
library(Biostrings)

taxonomy.file.path <- "/home/aubbxs/sh_general_release_dynamic_s_all_29.11.2022_Euk_SynMock.fasta"
otus.file.path <- "/scratch/aubbxs/fungi/21-22/output/cluster/otus_R1.fasta"

# Fasta
FASTA.otus <- readDNAStringSet(otus.file.path, format="fasta", seek.first.rec=TRUE, use.names=TRUE)

taxa <- assignTaxonomy(FASTA.otus, taxonomy.file.path, multithread=TRUE, tryRC = TRUE)

saveRDS(taxa, file = "taxa_out_DADA2_NBC.rds")
