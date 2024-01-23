# libraries
library(Biostrings)
library(dplyr)

# Load in otus.fasta
FASTA.fungi <- readDNAStringSet("otus_R1.fasta", seek.first.rec=TRUE, use.names=TRUE)
FASTA.fungi.dataframe <- as.data.frame(FASTA.fungi)
FASTA.fungi.dataframe$OTU <- rownames(FASTA.fungi.dataframe)

# Load in NBC taxonomy
NBC_taxonomy <- read.csv("NBCtaxa.csv")
head(NBC_taxonomy)

# NBC taxonomy does not contain the OTUs associated with our sequences so we use the otus.fasta and taxonomy to align OTU# and sequence together in a column called OTU_sequence
NBC_taxonomy2 <- left_join(FASTA.fungi.dataframe, NBC_taxonomy, by = c("x" = "X"))
head(NBC_taxonomy2)

#Edit the taxonomy file to remove unnecessary aspects. We clean the taxonomy table for ease of use.
nbc.tax.fungi <- as.data.frame(NBC_taxonomy2)

rownames(nbc.tax.fungi) <- nbc.tax.fungi$OTU
nbc.tax.fungi$Kingdom <- gsub('k__','', nbc.tax.fungi$Kingdom)
nbc.tax.fungi$Phylum <- gsub('p__','', nbc.tax.fungi$Phylum)
nbc.tax.fungi$Class <- gsub('c__','', nbc.tax.fungi$Class)
nbc.tax.fungi$Order <- gsub('o__','', nbc.tax.fungi$Order)
nbc.tax.fungi$Family <- gsub('f__','', nbc.tax.fungi$Family)
nbc.tax.fungi$Genus <- gsub('g__','', nbc.tax.fungi$Genus)
nbc.tax.fungi$Species <- gsub('s__','', nbc.tax.fungi$Species)

# Replace na with unidentified
nbc.tax.fungi <- replace(nbc.tax.fungi, is.na(nbc.tax.fungi), "unidentified")

# Move unidentified to a column called lowest taxonomic rank
nbc.tax.fungi$Lowest_Taxnomic_Rank <- ifelse(nbc.tax.fungi$Phylum == "unidentified", nbc.tax.fungi$Kingdom,
                                             ifelse(nbc.tax.fungi$Class == "unidentified", nbc.tax.fungi$Phylum,
                                                    ifelse(nbc.tax.fungi$Order == "unidentified", nbc.tax.fungi$Class,
                                                           ifelse(nbc.tax.fungi$Family == "unidentified", nbc.tax.fungi$Order,
                                                                  ifelse(nbc.tax.fungi$Genus == "unidentified", nbc.tax.fungi$Family,
                                                                         ifelse(nbc.tax.fungi$Species == "unidentified", nbc.tax.fungi$Genus, 
                                                                                paste(nbc.tax.fungi$Genus, nbc.tax.fungi$Species, sep = "_")))))))

nbc.tax.fungi$Label <- paste(nbc.tax.fungi$OTU, nbc.tax.fungi$Lowest_Taxnomic_Rank, sep = "_")

# Load in SINTAX taxonomy
SINTAXa <- read.csv("fungi_R1_UNITE2.csv", header = F, na.strings = "")
SINTAXa <- replace(SINTAXa, is.na(SINTAXa), "unidentified")
colnames(SINTAXa) = c("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
SINTAXa$Lowest_Taxonomic_Rank <- ifelse(SINTAXa$Phylum == "unidentified", SINTAXa$Kingdom,
                                             ifelse(SINTAXa$Class == "unidentified", SINTAXa$Phylum,
                                                    ifelse(SINTAXa$Order == "unidentified", SINTAXa$Class,
                                                           ifelse(SINTAXa$Family == "unidentified", SINTAXa$Order,
                                                                  ifelse(SINTAXa$Genus == "unidentified", SINTAXa$Family,
                                                                         ifelse(SINTAXa$Species == "unidentified", SINTAXa$Genus,
                                                                                paste(SINTAXa$Genus, SINTAXa$Species, sep = "_")))))))

SINTAXa$Label <- paste(SINTAXa$OTU, SINTAXa$Lowest_Taxonomic_Rank, sep = "_")

# Subset SINTAX to only the Mocki identified OTUs
mock.OTU <- SINTAXa[SINTAXa$Kingdom == "Mocki",]

# Left join to NBC taxonomy
nbc.tax.fungi0 <- left_join(nbc.tax.fungi, mock.OTU, by = "OTU")

OTU.mock <- mock.OTU$OTU

nbc.tax.fungi4 <- nbc.tax.fungi0

# Mocki will be subtituted in place of the NBC taxonomy that was previously unidentified (this will not happen for all samples, some are purely unidentified). If not a mock,
# it will remain the same as the original NBC taxonomy call
nbc.tax.fungi4$Kingdom.x <- ifelse(nbc.tax.fungi4$OTU %in% OTU.mock, nbc.tax.fungi4$Kingdom.y, nbc.tax.fungi4$Kingdom.x)
nbc.tax.fungi4$Phylum.x <- ifelse(nbc.tax.fungi4$OTU %in% OTU.mock, nbc.tax.fungi4$Phylum.y, nbc.tax.fungi4$Phylum.x)
nbc.tax.fungi4$Class.x <- ifelse(nbc.tax.fungi4$OTU %in% OTU.mock, nbc.tax.fungi4$Class.y, nbc.tax.fungi4$Class.x)
nbc.tax.fungi4$Order.x <- ifelse(nbc.tax.fungi4$OTU %in% OTU.mock, nbc.tax.fungi4$Order.y, nbc.tax.fungi4$Order.x)
nbc.tax.fungi4$Family.x <- ifelse(nbc.tax.fungi4$OTU %in% OTU.mock, nbc.tax.fungi4$Family.y, nbc.tax.fungi4$Family.x)
nbc.tax.fungi4$Genus.x <- ifelse(nbc.tax.fungi4$OTU %in% OTU.mock, nbc.tax.fungi4$Genus.y, nbc.tax.fungi4$Genus.x)
nbc.tax.fungi4$Species.x <- ifelse(nbc.tax.fungi4$OTU %in% OTU.mock, nbc.tax.fungi4$Species.y, nbc.tax.fungi4$Species.x)
nbc.tax.fungi4$Label.x <- ifelse(nbc.tax.fungi4$OTU %in% OTU.mock, nbc.tax.fungi4$Label.y, nbc.tax.fungi4$Label.x)

# Get rid of the SINTAX portion
nbc.tax.fungi5 <- nbc.tax.fungi4 %>%
  select(OTU:Label.x)
rownames(nbc.tax.fungi5) <- nbc.tax.fungi5$OTU

# Rename columns
colnames(nbc.tax.fungi5) = c("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Lowest_Taxonomic_Rank", "Label")

# Write csv
write.csv(nbc.tax.fungi5, file = "NBC_SINTAX_comb.csv")
