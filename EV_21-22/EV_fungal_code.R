library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(ggpubr)
library(decontam)
library(metagenomeSeq)
library(indicspecies)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Taxonomy 
#Edit the taxonomy rds
nbc.tax.fungi <- readRDS("./21-22/taxa_out_DADA2_NBC.rds")

nbc.tax.fungi <- as.data.frame(nbc.tax.fungi)
nbc.tax.fungi$OTU <- sprintf("FOTU_%s", seq(1:nrow(nbc.tax.fungi)))
rownames(nbc.tax.fungi) <- nbc.tax.fungi$OTU
nbc.tax.fungi$Kingdom <- gsub('k__','', nbc.tax.fungi$Kingdom)
nbc.tax.fungi$Phylum <- gsub('p__','', nbc.tax.fungi$Phylum)
nbc.tax.fungi$Class <- gsub('c__','', nbc.tax.fungi$Class)
nbc.tax.fungi$Order <- gsub('o__','', nbc.tax.fungi$Order)
nbc.tax.fungi$Family <- gsub('f__','', nbc.tax.fungi$Family)
nbc.tax.fungi$Genus <- gsub('g__','', nbc.tax.fungi$Genus)
nbc.tax.fungi$Species <- gsub('s__','', nbc.tax.fungi$Species)

nbc.tax.fungi <- replace(nbc.tax.fungi, is.na(nbc.tax.fungi), "unidentified")

nbc.tax.fungi$Lowest_Taxnomic_Rank <- ifelse(nbc.tax.fungi$Phylum == "unidentified", nbc.tax.fungi$Kingdom,
                                             ifelse(nbc.tax.fungi$Class == "unidentified", nbc.tax.fungi$Phylum,
                                                    ifelse(nbc.tax.fungi$Order == "unidentified", nbc.tax.fungi$Class,
                                                           ifelse(nbc.tax.fungi$Family == "unidentified", nbc.tax.fungi$Order,
                                                                  ifelse(nbc.tax.fungi$Genus == "unidentified", nbc.tax.fungi$Family,
                                                                         ifelse(nbc.tax.fungi$Species == "unidentified", nbc.tax.fungi$Genus, 
                                                                                paste(nbc.tax.fungi$Genus, nbc.tax.fungi$Species, sep = "_")))))))

nbc.tax.fungi$Label <- paste(nbc.tax.fungi$OTU, nbc.tax.fungi$Lowest_Taxnomic_Rank, sep = "_")

#create taxonomy variable
tax <- nbc.tax.fungi
head(tax)
tax$OTU <- rownames(tax)
TAX.fungi <- phyloseq::tax_table(as.matrix(tax))

# OTU Table 
table <- read.csv("./21-22/otu.table.csv")
rownames(table) <- table$OTU
table <- table[,-1]
OTU.fungi <- phyloseq::otu_table(table, taxa_are_rows = TRUE)

# Metadata 
samples <- read.csv("./21-22/21-22_Metadata.csv", na.strings = "na")
rownames(samples) <- samples$Identifier #row names must match OTU table headers
SAMP.fungi <- phyloseq::sample_data(samples)

# Fasta File
FASTA.fungi <- Biostrings::readDNAStringSet("./21-22/otus_R1.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

# Phyloseq Object Creation
phyloseq.start <- phyloseq(SAMP.fungi, TAX.fungi, OTU.fungi, FASTA.fungi)

#save phyloseq object as an rds
saveRDS(phyloseq.start, "21-22-fungi-phyloseq.rds")

# New number of total reads
sum(sample_sums(EVSamples))

# Mean and median read depth 
mean(sample_sums(EVSamples))

median(sample_sums(EVSamples))

# Histogram including median read depth
read.depths <- data.frame(sample_sums(EVSamples))
colnames(read.depths) <- "read.depth"
read.depth.plot <- ggplot(read.depths, aes(read.depth)) +
  geom_histogram(fill = cbbPalette[[3]], color = "black") + 
  geom_vline(xintercept = median(sample_sums(EVSamples)), linetype = "dashed") + 
  theme_classic() + 
  xlab("Read Depth")
read.depth.plot

#Rarefaction

sam.data <- data.frame(EVSamples@sam_data)
sam.data$Sample_ID <- sam.data$Identifier
pOTU.table <- EVSamples@otu_table
S <- specnumber(t(pOTU.table)) # observed number of species
raremax <- min(rowSums(t(pOTU.table)))
Srare <- rarefy(t(pOTU.table), raremax)
rare.fun <- rarecurve(t(pOTU.table), step = 1000, sample = raremax, col = "blue", cex = 0.6)

prok.rare.curve.extract <- NULL
for(i in 1:length(rare.fun)){
  sample.200 <- data.frame(rare.spec = rare.fun[[i]])
  sample.200$read_depth <- attr(rare.fun[[i]], "Subsample")
  sample.200$Sample_ID <- rownames(t(pOTU.table[,i]))
  prok.rare.curve.extract <- rbind.data.frame(prok.rare.curve.extract, sample.200)
}
prok.rare.curve.extract2 <- left_join(sam.data, prok.rare.curve.extract, by = "Sample_ID")

rare.curve <- ggplot(prok.rare.curve.extract2, aes(x = read_depth, y = rare.spec, group = Sample_ID)) + 
  #geom_point() +
  geom_line(color = "grey") + 
  xlab("Reads") + 
  ylab("Number of OTUs") + 
  theme_classic() + 
  geom_vline(xintercept = median(sample_sums(EVSamples)), linetype = "dashed")
rare.curve

#Alpha Diversity

EVSamples@sam_data$shannon <- estimate_richness(EVSamples, measures=c("Shannon"))$Shannon
EVSamples@sam_data$invsimpson <- estimate_richness(EVSamples, measures=c("InvSimpson"))$InvSimpson
EVSamples@sam_data$richness <- estimate_richness(EVSamples, measures=c("Observed"))$Observed
EVSamples@sam_data$even <- EVSamples@sam_data$shannon/log(EVSamples@sam_data$richness)

sample.data.fungi <- data.frame(EVSamples@sam_data)

#CSS Normalization

MGS <- phyloseq_to_metagenomeSeq(EVSamples)
p <- metagenomeSeq::cumNormStatFast(MGS)
MGS <- metagenomeSeq::cumNorm(MGS, p =p)
metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample
norm.fungi <- metagenomeSeq::MRcounts(MGS, norm = T)
norm.fungi.OTU <- phyloseq::otu_table(norm.fungi, taxa_are_rows = TRUE)

physeq.css <- phyloseq::phyloseq(norm.fungi.OTU, SAMP.fungi, TAX.fungi, FASTA.fungi)

#Beta Diversity (PCoA)

#Subset
EVCSS <- subset_samples(physeq.css, EV.Time.Point %in% c("T6"))

# Principle coordinates analysis with Bray-Curtis distances
ordination.pcoa <- ordinate(EVCSS, "PCoA", "bray") # calculate the resemblance and ordinate using PCoA
ordination.pcoa$vectors # positions of your points on the PCoA graph
ordination.pcoa$values #values to calculate the variance explained on each axis (dimension)

pcoa <- plot_ordination(EVCSS, ordination = ordination.pcoa, type = "samples", color = "Fungicide") +
  theme_classic() + 
  geom_point(size = 4) +
  scale_color_manual(values = cbbPalette)
pcoa

#PERMANOVA

prok.dist.bray = phyloseq::distance(EVCSS, "bray") # create bray-curtis distance matrix
adonis2(prok.dist.bray~Fungicide*Cultivar, as(sample_data(EVCSS), "data.frame"))

#Dispersion

EVdisp <- betadisper(prok.dist.bray, EVCSS@sam_data$Fungicide)
permutest(EVdisp)

#convert rds to csv

write.csv(nbc.tax.fungi, "NBC.csv")

