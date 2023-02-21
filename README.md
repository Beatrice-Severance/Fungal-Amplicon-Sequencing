# Fungal-Amplicon-Sequencing

A pipeline that encompasses the computational steps of ITS amplicon sequencing, starting from demultiplexed reads.
The scripts in this pipeline were originally run on the [Alabama supercomputer](https://www.asc.edu/) and thus might make reproducibility difficult for users who do not utilize this system.
Scripts for this pipeline are grouped in a folder along with the "2022samples.txt" file. The Alabama HPC utilizes a slurm queue system where jobs are submitted and run. The scripts being located in the same file along with the samples file helps run jobs in the queue more efficiently.

# Strip Primers
This is a loop script that will remove primers that may be included in the output from the merging reads step. The "2022samples.txt" file will be used in this step as well. Only forward (R1) reads are used in the fungal pipeline because ITS1 is the region of interest. The primer for this region is:
- CTTGGTCATTTAGAGGAAGTAA

# Run Statistics
This is a script that uses VSEARCH to run statistics on their dataset.

# Filter and Trim
This is a script that will filter out bad quality reads from the previous step. Parameters are set at an e-value of 1.0 and a length of 263bp. Parameters can be edited based on user needs. The left side of the sequences are trimmed for the purposes of the analysis, and FastQC is run on the data to provide statistics. Users can determine whether they need to filter/trim further based on the results.

# Cluster
This is a script that can dereplicate, cluster, and remove chimeras.
- De-noising step will provide zero-radius OTUs (zOTUs).
- Clustering will provide OTUs based on traditional 97% identity.
- USEARCH is a program that is utilized for the de-noising and clustering steps. For more information on these programs the following links can be used:
- [UPARSE vs. UNOISE](http://www.drive5.com/usearch/manual/faq_uparse_or_unoise.html)
- [otutab command](http://www.drive5.com/usearch/manual/cmd_otutab.html)
- [Sample identifiers in read labels](http://www.drive5.com/usearch/manual/upp_labels_sample.html)
- [Bugs and fixes for USEARCH v11](http://drive5.com/usearch/manual/bugs.html)
- [Technical support](http://drive5.com/usearch/manual/support.html) 

# Create Taxonomy
This script utilizes the NBC algorithm in order to produce a fungal taxonomy from the provided samples. This script requires the file dada2_assigntax_NBC.R which will use the database provided to create NBC taxonomy. This script requires a large amount of memory in order to run, likely because the database itself is large. The script will create a .rds file which will allow users to view the taxonomy in R or RStudio.

# Mapping
This script will create an OTU table that will be used for downstream analysis. It utilizes the demultiplexed reads and aligns these reads back to the clustered OTUs.

#
Combined, these steps should provide output files that can be utilized in a phyloseq object in R.

# R Analysis
