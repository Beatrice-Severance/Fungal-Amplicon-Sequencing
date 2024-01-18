# Fungal-Amplicon-Sequencing

A pipeline that encompasses the computational steps of ITS amplicon sequencing, starting from demultiplexed reads.
The scripts in this pipeline were originally run on the [Alabama supercomputer](https://www.asc.edu/) and thus might make reproducibility difficult for users who do not utilize this system. The Alabama supercomputer recently shifted from DMC to ASAX. Currently, users can utilize the [ASAX scripts](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/tree/main/ASAX_Scripts) to perform data processing, and these will be linked throughout the README file below. For reference, the DMC scripts previously used are linked [here](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/tree/main/DMC_Scripts).
Scripts for this pipeline are grouped in a folder along with the sample files from [2021](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/HPC_Scripts/2021_EV_fungal_samples.txt) and [2022](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/HPC_Scripts/2022_EV_fungal_samples.txt) leaf sequences extracted from pecan trees at the E.V. Smith facility in Alabama. The Alabama HPC utilizes a slurm queue system where jobs are submitted and run. The scripts being located in the same file along with the sample files help run jobs in the queue more efficiently.

# Creating a [samples.txt](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/DMC_Scripts/combsamples.txt) File
In order for loops to be performed throughout the pipeline, a file was created to help the computer quickly identify files and how to rename them. This shortens sequence files from looking like this:
```
132BFT6_S134_L001_R1_001.fastq.gz
112CCT1_S186_L001_R1_001.fastq.gz
112CCT2_S157_L001_R1_001.fastq.gz
112CCT3_S266_L001_R1_001.fastq.gz
```
to this
```
132BFT6_S134
112CCT1_S186
112CCT2_S157
112CCT3_S266
```
This samples file will be located in the scripts folder, since it makes it easy for the computer to find its reference quickly.

# [Strip Primers](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/ASAX_Scripts/1_strippingPrimers.sh)
This is a loop script that will remove primers that may be included in the demultiplexed reads. The sample files will be used in this step. Only forward (R1) reads are used in this fungal pipeline because ITS1 is the region of interest. The primer for this region is:
- CTTGGTCATTTAGAGGAAGTAA

# Combine Data
Before running statistics and filtering, it is necessary to combine the trimmed data from both the previous step. Stripped FASTQ files (*_trimmed.fastq) are concatenated into one file, trimmed_combined.fastq.

# [Run Statistics](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/ASAX_Scripts/2_stat_results.sh)
This is a script that uses VSEARCH to run statistics on the trimmed_combined.fastq file. Stats can be viewed to determine filtering parameters.

# [Filter](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/ASAX_Scripts/3_filtering_and_trimming.sh)
This is a script that utilizes VSEARCH to filter out bad quality reads from the trimmed_combined.fastq file. Parameters are set at an e-value of 1.0 and a length of 263bp, parameters determined from the previous statistics step. The left side of the sequences are trimmed by 44 bases with the -fastq_stripleft command. FastQC is run on the data to provide statistics. Users can determine whether they need to filter further based on the results.

# [Cluster](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/ASAX_Scripts/4_clustering.sh)
This is a script that utilizes VSEARCH to dereplicate, and USEARCH to cluster and remove chimeras.
- De-noising step will provide zero-radius OTUs (zOTUs).
- Clustering will provide OTUs based on traditional 97% identity.
- USEARCH is a program that is utilized for the de-noising and clustering steps. For more information on these programs the following links can be used:
- [UPARSE vs. UNOISE](http://www.drive5.com/usearch/manual/faq_uparse_or_unoise.html)
- [otutab command](http://www.drive5.com/usearch/manual/cmd_otutab.html)
- [Sample identifiers in read labels](http://www.drive5.com/usearch/manual/upp_labels_sample.html)
- [Bugs and fixes for USEARCH v11](http://drive5.com/usearch/manual/bugs.html)
- [Technical support](http://drive5.com/usearch/manual/support.html) 

# Create Taxonomy
Fungal taxonomy is created utilizing the UNITE database. Multiple taxonomy programs are used, SINTAX and NBC (Naive Baysian Classifier). SINTAX is able to classify mock sequences used in the experiment that the NBC commonly misses. However, SINTAX does not classify some well-known sequences, so the NBC algorithm is used to fill in the gaps. Both programs use the otus_R1.fasta file created in the clustering step as input. The NBC file was used as the base for eventual R integration in a phyloseq object after manually adding SINTAX mock sequences into corresponding unidentified NBC rows.
The mock sequences used come from the following [manuscript](https://doi.org/10.7717%2Fpeerj.4925).

## [NBC Taxonomy](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/ASAX_Scripts/5_DADA2_NBC.sh)
This script requires the R script [dada2_assigntax_NBC.R](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/ASAX_Scripts/dada2_assigntax_NBC.R) which will use the database provided to create taxonomy based on the NBC algorithm. This script requires a large amount of memory in order to run, so should be taken into account when attempting to replicate this pipeline. 1 core and 32gb of memory produced parallel efficiency of 99.37% and memory efficiency of 60.63%. The script creates an .rds file which will allow users to view the taxonomy in R or RStudio. The database that was used for this script is located [here](https://doi.plutof.ut.ee/doi/10.15156/BIO/2483914). The following release was used for analysis: 29.11.2022
Mock sequences in NBC format can be found [here](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Mock_Sequences/mocksequencesNBC.txt).

## [SINTAX Taxonomy](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/ASAX_Scripts/5_taxonomy_SINTAX.sh)
This script utilizes the SINTAX algorithm to create a fungal taxonomy. VSEARCH is the medium used to achieve this goal. The database that was used for this script is located [here](https://doi.plutof.ut.ee/doi/10.15156/BIO/2483924). The following release was used for analysis: 29.11.2022
Mock sequences in SINTAX format can be found [here](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Mock_Sequences/mocksequencesSINTAX.txt).

# [Mapping](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/ASAX_Scripts/6_mapping.sh)
This script will create an OTU table that will be used for downstream analysis. It utilizes the [demultiplexed reads](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Scripts/combsamples.txt) and aligns these reads back to the clustered OTUs (otus_R1.fasta). Before the main part of the script is run, a [Python script](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/ASAX_Scripts/replacefastaheaders_filename.py) is used to replace fasta headers. This is performed in order for the OTUs to be mapped back to their respective samples, instead of being processed as a singular sample. 

#
Combined, these steps provide the following output files that can be utilized in a phyloseq object in R:
- [OTU table](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/EV_21-22/EV_21-22_phyloseq_input/otu.table.csv)
- ITS taxonomy files (from SINTAX and [NBC](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/EV_21-22/EV_21-22_phyloseq_input/NBC.csv))
- [otus.fasta](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/EV_21-22/EV_21-22_phyloseq_input/otus_R1.fasta) file

# R Analysis
R Analysis for this project starts with the creation of a phyloseq object. In addition to the files above, a [metadata](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/EV_21-22/EV_21-22_phyloseq_input/21-22_Metadata.csv) file will be necessary to run analysis. The [R Markdown file](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/EV_21-22/EV_21-22_Fungi.Rmd) will execute the following steps:
- Load dependencies
  - Dependencies used for analysis:
    - phyloseq [version 1.44.0](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217#s6) 
    - vegan [version 2.6-4](https://github.com/vegandevs/vegan/releases/tag/v2.6-4)
    - tidyverse [version 2.0.0](https://github.com/tidyverse/tidyverse/releases/tag/v2.0.0)
    - ggplot2 [version 3.4.2](https://cloud.r-project.org/web/packages/ggplot2/index.html)
    - Biostrings [version 2.68.1](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
    - ggpubr [version 0.6.0](https://cran.r-project.org/web/packages/ggpubr/index.html)
    - decontam [version 1.20.0](https://github.com/benjjneb/decontam)
    - metagenomeSeq [version 1.42.0](https://github.com/HCBravoLab/metagenomeSeq)
    - indicspecies [version 1.7.14](https://cran.r-project.org/web/packages/indicspecies/index.html)
- Utilize a colorblind palette
- Load in [files](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/tree/main/EV_21-22/EV_21-22_phyloseq_input) to create a phyloseq object
- Decontaminate samples 
  - Take out controls
  - Remove low-quality reads (<5000)
  - Subset to kingdom Fungi
- Provide read distribution for the dataset, including a histogram
- Rarefaction analysis, including line graphs
- Alpha diversity analysis
  - Shannon
  - Inverse Simpson
  - Richness
- Cumulative sum scaling (CSS) Normalization
- Beta diversity analysis
  - Principal coordinates analysis (PCoA) plot with Bray-Curtis distances
- PERMANOVA to test for differences in centroids

The rendered HTML file for the data can be viewed [here](https://htmlpreview.github.io/?https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/EV_21-22_Fungi.html).

Figures generated in R can be viewed [here](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/tree/main/EV_21-22/EV_21-22_figures).