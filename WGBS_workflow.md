# WGBS Workflow Overview

This document outlines the steps performed by the Snakemake workflow in `workflow/snakefile`. The input files must be 
whole genome bisulfite sequence files in FASTQ format. The final output files are filtered SNV files. 

This workflow is heavily based off of the pipelines in the [Github](https://github.com/MLindner0/lindner_et_al-2021-mer-snps_from_bs_data) repository for [Lindner *et al.,* 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13493). 

## Sample Pre-processing

### 1. Split files for faster processing
**Rule split**

Raw fastq files are split into user-defined chunks using the partition.sh 
script as part of BBMap. Files are normally split so that resulting files 
contain fewer than 100 million reads.

### 2. Quality trimming with TrimGalore
**Rule trimgalore**

Reads are trimmed in each split file using the program TrimGalore, which is 
a wrapper script for Cutadapt. 

## Read alignment
### 3. Prepare genome
**Rule prep_genome**

The reference genome is prepared for alignment with the Bismark tool 
`bismark_genome_preparation`. By default, Bismark aligns using Bowtie2, so this 
step creates a Bowtie2-indexed reference genome. 

### 4. Align reads
**Rule align**

Bismark is used to align the reads to the reference genome. Because the SNV caller 
CGMapTools is used downstream, Bismark is run with the `--old_flag` option. 

### 5. Merge split bam files
**Rule merge_bams**

Bam files generated from part of the trimmed fastq files are merged using 
Picard MergeSamFiles. Note that multiplexed samples are still separate.

### 6. Sequence deduplication
**Rule dedup**

The Bismark tool `deduplicate_bismark` is used to remove duplicate reads that can 
arise from PCR amplification.

### Add read group information
**Rule add_rgids**

Picard AddOrReplaceReadGroups is used to add read group information to each file. 
Read group information is derived from `resources/reference_files/rg_labels.txt`. 

### 7. Merge aligned files by sample
**Rule merge**

Picard MergeSamFiles is used to merge aligned SAM files from sample libraries 
sequenced on different flow cell lanes.

### 8. Sort merged files
**Rule sortsam**

Picard SortSam is used to sort the merged aligned files by coordinate. Files were 
previously name-sorted as this is required for Bismark deduplication.

### 9. Index sorted files
**Rule index**

BAM file indexes are created using samtools.

## Statistics
### 10. Calculate depth and breadth of read coverage
**Rule coverage**

Samtools depth is used to calculate the depth of coverage, samtools mpileup is used 
to calculate the breadth of coverage (>4x) across the genome, and samtools flagstat 
returns read mapping statistics. 

## SNV calling
### 11. CGMapTools conversion
**Rule cgmap_conversion**

Prior to SNV calling, all aligned BAM files must be transformed into ATCGmap files. 
This is done with the `CGMapTools convert bam2cgmap` utility.

### 12. Variants are called on each sample
**Rule cgmap_bayes**

SNPs are called using CGMapTools in high-precision Bayesian mode (`-m bayes --bayes-dynamicP`).


