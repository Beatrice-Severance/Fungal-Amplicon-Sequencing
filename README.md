# Fungal-Amplicon-Sequencing

A pipeline that encompasses the computational steps of ITS amplicon sequencing, starting from demultiplexed reads.
The scripts in this pipeline were originally run on the [Alabama supercomputer](https://www.asc.edu/) and thus might make reproducibility difficult for users who do not utilize this system.
Scripts for this pipeline are grouped in a folder along with the sample files from [2021](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Scripts/2021_EV_fungal_samples.txt) and [2022](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Scripts/2022_EV_fungal_samples.txt) leaf sequences extracted from pecan trees at the E.V. Smith facility in Alabama. The Alabama HPC utilizes a slurm queue system where jobs are submitted and run. The scripts being located in the same file along with the sample files help run jobs in the queue more efficiently.

# [Strip Primers](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Scripts/stripping_primers.sh)
This is a loop script that will remove primers that may be included in the demultiplexed reads. The sample files will be used in this step. Only forward (R1) reads are used in this fungal pipeline because ITS1 is the region of interest. The primer for this region is:
- CTTGGTCATTTAGAGGAAGTAA

# Combine 2021 and 2022 Data
Before running statistics and filtering, it is necessary to combine data from both 2021 and 2022. Stripped FASTQ files (*_trimmed.fastq) are concatenated into one file, trimmed_combined.fastq.

# [Run Statistics](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Scripts/stats.sh)
This is a script that uses VSEARCH to run statistics on the trimmed_combined.fastq file. Stats can be viewed to determine filtering parameters.

# [Filter](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Scripts/filter_and_trim.sh)
This is a script that utilizes VSEARCH to filter out bad quality reads from the trimmed_combined.fastq file. Parameters are set at an e-value of 1.0 and a length of 263bp, parameters determined from the previous statistics step. The left side of the sequences are trimmed by 44 bases with the -fastq_stripleft command. FastQC is run on the data to provide statistics. Users can determine whether they need to filter further based on the results.

# [Cluster](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Scripts/clustering.sh)
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
Fungal taxonomy is created utilizing the UNITE database. Multiple taxonomy programs are used, SINTAX and NBC (Naive Baysian Classifier). SINTAX is able to classify mock sequences used in the experiment that the NBC commonly misses. However, SINTAX does not classify some well-known sequences, so the NBC algorithm is used to fill in the gaps. Both programs use the otus_R1.fasta file created in the clustering step as input.
The mock sequences used come from the following [manuscript](https://doi.org/10.7717%2Fpeerj.4925).

## [NBC Taxonomy](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Scripts/5_DADA2_NBC.sh)
This script requires the file dada2_assigntax_NBC.R which will use the database provided to create taxonomy based on the NBC algorithm. This script requires a large amount of memory in order to run, so should be taken into account when attempting to replicate this pipeline. 1 core and 32gb of memory produced parallel efficiency of 99.37% and memory efficiency of 60.63%. The script creates an .rds file which will allow users to view the taxonomy in R or RStudio. The database that was used for this script is located [here](https://doi.plutof.ut.ee/doi/10.15156/BIO/2483914). The latest release was used for analysis: 29.11.2022
Mock sequences in NBC format can be found [here](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Mock_Sequences/mocksequencesNBC.txt).

## [SINTAX Taxonomy](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Scripts/5_taxonomy_SINTAX.sh)
This script utilizes the SINTAX algorithm to create a fungal taxonomy. VSEARCH is the medium used to achieve this goal. The database that was used for this script is located [here](https://doi.plutof.ut.ee/doi/10.15156/BIO/2483924). The latest release was used for analysis: 29.11.2022
Mock sequences in SINTAX format can be found [here](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Mock_Sequences/mocksequencesSINTAX.txt).

# [Mapping](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Scripts/6_mapping.sh)
This script will create an OTU table that will be used for downstream analysis. It utilizes the [demultiplexed reads](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/Scripts/combsamples.txt) and aligns these reads back to the clustered OTUs (otus_R1.fasta).

#
Combined, these steps provide the following output files that can be utilized in a phyloseq object in R:
- [OTU table](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/phyloseq_input/otu.table.csv)
- ITS taxonomy files (from SINTAX and [NBC](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/phyloseq_input/NBC.csv))
- [otus.fasta](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/phyloseq_input/otus_R1.fasta) file

# R Analysis
R Analysis for this project starts with the creation of a phyloseq object. In addition to the files above, a [metadata](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/phyloseq_input/21-22_Metadata.csv) file will be necessary to run analysis. The [R Markdown file](https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/EV_21-22_Fungi.Rmd) will execute the following steps:
- Load dependencies
- Utilize a colorblind palette
- Load above files to create a phyloseq object
- Decontaminate samples (take out controls, low quality reads, etc.)
- Provide read distribution for the dataset, including a histogram
- Rarefaction analysis, including line graphs
- Alpha diversity analysis
- Cumulative sum scaling (CSS) Normalization
- Beta diversity analysis, including a principal coordinates analysis (PCoA) plot with Bray-Curtis distances
- PERMANOVA to test for differences in centroids

The rendered HTML file for the data can be viewed [here](https://htmlpreview.github.io/?https://github.com/Beatrice-Severance/Fungal-Amplicon-Sequencing/blob/main/EV_21-22_Fungi.html).