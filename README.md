# variant-calling

author: Sam Lee

code provided as-is where-is etc

Simple variant calling pipeline for RNA-seq data. Implemented using Opossum and Platypus for data preparation and variant calling.

A **VCF** of all variants and a  gzipped **VCF** per chromosme will be produced.

## Usage

```
# install Snakemake if needed 
#  (easiest through conda)
#  (all other dependencies managed by conda)
conda create -n myEnv Snakemake
source activate myEnv

Snakemake --use-conda --cores N

```

In `config.yaml` **bamDir** should point to the directory with the aligned + sorted bam files for the analysis. 

**samples** is then a dictionary with *key*:*value* pairs of *sample*:*sample-bam-file* although the value doesn't really matter as the Snakefile expects the input to be of the format `{sample}_Aligned.sortedByCoord.out.bam`

## To do

* Write rule for mapping of fastq files with STAR
* Decide on imputation method and implement
