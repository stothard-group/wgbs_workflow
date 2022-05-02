# Whole Genome Bisulfite Sequencing Workflow

 - FILE: workflow/snakefile
 - AUTH: Emily Herman (eherman@ualberta.ca)
 - DATE: MAY 5, 2022
 - VERS: 1.0

This Snakemake workflow processes whole genome bisulfite sequence data and extracts single nucleotide variants. The 
input is short read WGBS data in fastq format, and the output is a genotyped VCF file for each sample.

README Sections

[Installation](#installation)

[Setup](#setup)

[Input Files](#input-files)

[Running the Workflow](#running-the-workflow)

[Memory Requirements](#memory-requirements)

[Output Files](#output-files)



## Installation

Download the repository

```angular2html
git clone https://github.com/stothard-group/wgbs_workflow.git
```

## Setup

Requirements:

 - [Snakemake](https://snakemake.readthedocs.io/en/stable/)
 - [Python3](https://www.python.org/downloads/)
 - [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
 - [CGMapTools](https://cgmaptools.github.io/) - Note that CGMapTools requires Python 2

The absolute paths to the programs TrimGalore and CGMapTools must be supplied in the `config/WGBS.config.yaml` file.

If the workflow is not run on a Compute Canada cluster, install the following programs and remove `module load ...` 
lines from shell commands in the `workflow/snakefile` file. Programs must be installed in the user's path, or their 
path must be given in the shell command. Note that many of these tools require their own dependencies.

 - [BBMap](https://github.com/BioInfoTools/BBMap) (specifically partition.sh)
 - [Java SE Development Kit>=14.0.2](https://www.oracle.com/ca-en/java/technologies/javase/jdk13-archive-downloads.html)
 - [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
 - [Picard=2.26.3](https://github.com/broadinstitute/picard)
 - [Samtools](http://www.htslib.org/)
 - [Genome Analysis Toolkit>=4.2.4.0](https://github.com/broadinstitute/gatk)
 - [Bcftools](https://samtools.github.io/bcftools/)


Install the following python modules:
 - pandas >= 1.3.5

## Input files

The following inputs are required:
 - Directory containing raw WGBS fastq files
 - A reference genome
 - A tab-formatted file containing read group information with the fields ID, SM, LB, PL for all samples. [See this 
GATK article on read groups for more information](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)
 - A file containing a list of contig names; one name per line

Raw fastq files should be gzip-compressed with the filename structure `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`. 
Read group ID information will be used to merge multiplexed samples.

The paths to the required input files and directories should be edited in the file `config/WGBS.config.yaml`. This file 
also contains species- and analysis-specific parameters affecting workflow function.



---
*NOTE*

If using NovaSeq6000 samples, add `--2colour 20` to the shell command for `rule trimgalore`.

---


## Running the workflow

The following command will run up to 20 workflow jobs on a SLURM cluster:

snakemake --snakefile workflow/snakefile --cluster "sbatch -N 1 -c {resources.cores} --mem={resources.mem_mb} 
--time={resources.runtime} --account=account" --rerun-incomplete --printshellcmds -j 20


## Memory requirements

This workflow requires at most 16 cores and 64GB RAM for running Bismark alignment and deduplication steps. The longest 
job is Bismark alignment, for which 12 hours is allocated by default. Resource limits can be adjusted for individual 
rules within the `workflow/snakefile` file.


---
*NOTE*

The core and RAM requirements for Bismark alignment (`rule align`) have a special relationship to the value of the 
`threads` parameter. Please consult the information on the `--parallel` option in [Bismark documentation Appendix (II)](https://github.com/FelixKrueger/Bismark/tree/master/Docs) 
before changing any of these values.

---

## Output files

All intermediate and final results files are written to the `results/` directory. For each sample, a quality-filtered 
genotyped VCF file in a GATK-like format is produced: `snps/{sample}.goodSNPs.vcf`.

Below is a list of all output directories and files created by the workflow in `results/`. Note that for non-multiplexed 
samples, {inputID} and {sample} values will be identical. Otherwise, {sample} name is determined by the read group ID 
file specified in `config/WGBS.config.yaml`.


| Directory          | File name                                     | Description                                                 |
|--------------------|-----------------------------------------------|-------------------------------------------------------------|
| split              | {inputID}_R{1,2}_part{n}.fastq.gz             | Raw fastq files split into n-parts <br/>(specified in config)    |
| trimmed            | {inputID}_R{1,2}_part{n}_val_1.fq.gz          | Trimmed fastq files split into n-parts                      |
| align              | {inputID}_R1_part{n}_val_1_bismark_bt2_pe.bam | Bam alignment file produced by Bismark <br/>for each n-part      |
| align              | {inputID}.bismark.m.bam                       | Merged bam file for each input ID                           |
| base_quality_recal | {inputID}.deduplicated.bam                    | Bam file following read deduplication                       |
| base_quality_recal | {inputID}.deduplicated.rgid.bam               | Deduplicated bam file with read group IDs added             |
| base_quality_recal | {sample}.deduplicated.rgid.m.bam              | Bam files merged per sample if multiplexed                  |
| base_quality_recal | {sample}.deduplicated.rgid.m.sorted.bam       | Sorted merged bam file                                      |
| stat_and_coverage  | {sample}.average_coverage.txt                 | Average sequencing depth across the sample                  |
| stat_and_coverage  | {sample}.bases_covered.txt                    | Number of bases above the minimum depth <br/>specified in config |
| stat_and_coverage  | {sample}.alignment_stats.txt                  | Read alignment statistics                                   |
| cgmap              | {sample}.ATCGmap.gz                           | ATCGmap file derived from sample bam                        |
| cgmap              | {sample}.bayes_dynamicP.SNPs.{vcf,snv}        | Variants in both VCF and SNV format files                   |
| cgmap              | {sample}.bayes_dynamicP.SNPs.reheader.vcf     | VCF file with correct sample info in header                 |
| snps               | {sample}.bayes_dynamicP.SNPs.rh.noNs.vcf      | VCF file with non-A/T/G/C variants removed                  |
| snps               | {sample}.rawSNPs.vcf                          | VCF file containing only biallelic SNPs                     |
| snps               | {sample}.filtSNPs.vcf                         | Soft-filtered VCF file with FAIL_DP filter tag              |
| snps               | {sample}.goodSNPs.vcf                         | VCF file containing quality-filtered SNPs                   |


Benchmark files are created for every rule with cluster submission detailing job resource use. These files are found in 
the directory `benchmark/` and have the structure `{rule}.{inputFileBasename}.txt`.