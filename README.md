# Fungal-Amplicon-Sequencing

A pipeline that encompasses the computational steps of ITS amplicon sequencing, starting from demultiplexed reads.
The scripts in this pipeline were originally run on the [Alabama supercomputer](https://www.asc.edu/) and thus might make reproducibility difficult for users who do not utilize this system.
Scripts for this pipeline are grouped in a folder along with the sample files from 2021 and 2022 leaf sequences extracted from pecan trees at the E.V. Smith facility in Alabama (2021_EV_fungal_samples.txt and 2022_EV_fungal_samples.txt). The Alabama HPC utilizes a slurm queue system where jobs are submitted and run. The scripts being located in the same file along with the sample files help run jobs in the queue more efficiently.

# Strip Primers
This is a loop script that will remove primers that may be included in the demultiplexed reads. The sample files will be used in this step. Only forward (R1) reads are used in this fungal pipeline because ITS1 is the region of interest. The primer for this region is:
- CTTGGTCATTTAGAGGAAGTAA

# Combine 2021 and 2022 Data
Before running statistics and filtering, it is necessary to combine data from both 2021 and 2022. Stripped FASTQ files (*_trimmed.fastq) are concatenated into one file, trimmed_combined.fastq.

# Run Statistics
This is a script that uses VSEARCH to run statistics on the trimmed_combined.fastq file. Stats can be viewed to determine filtering parameters.

# Filter
This is a script that utilizes VSEARCH to filter out bad quality reads from the trimmed_combined.fastq file. Parameters are set at an e-value of 1.0 and a length of 263bp, parameters determined from the previous statistics step. The left side of the sequences are trimmed by 44 bases with the -fastq_stripleft command. FastQC is run on the data to provide statistics. Users can determine whether they need to filter further based on the results.

# Cluster
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

## NBC Taxonomy
This script requires the file dada2_assigntax_NBC.R which will use the database provided to create taxonomy based on the NBC algorithm. This script requires a large amount of memory in order to run, so should be taken into account when attempting to replicate this pipeline. 1 core and 32gb of memory produced parallel efficiency of 99.37% and memory efficiency of 60.63%. The script creates an .rds file which will allow users to view the taxonomy in R or RStudio.

## SINTAX Taxonomy
This script utilizes the SINTAX algorithm to create a fungal taxonomy. VSEARCH is the medium used to achieve this goal. 

# Mapping
This script will create an OTU table that will be used for downstream analysis. It utilizes the demultiplexed reads (combsamples.txt) and aligns these reads back to the clustered OTUs (otus_R1.fasta).

#
Combined, these steps provide output files that can be utilized in a phyloseq object in R.

# R Analysis
