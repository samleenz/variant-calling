
# 2018-03-08
# Sam Lee

configfile: "config.yaml"
from os.path import join
import re
#

chrNs = ["1","2","3","4","5","6","7","8","9","10","11","12","13",
           "14","15","16","17","18","19","20","21","22","X","Y"]

bam_list = expand("out/bam/{sample}_Aligned.sortedByCoord.out.Opossum.bam", sample=config["samples"])

rule all:
  input:
    expand("out/chr_compress/chr_{chr}_variants.sorted.vcf.gz", chr=chrNs)

# rule zip:
#   input:
#     "out/chromosomes/chr_{chr}_variants.sorted.vcf"
#   output:
#     "out/chr_compress/chr_{chr}_variants.sorted.vcf.gz"
#   shell:
#     """
#     cat {input} | bgzip -c > {output}
#     """

rule sort_vcf:
  input:
    "out/chromosomes/chr_{chr}_variants.recode.vcf"
  output:
    # "out/chromosomes/chr_{chr}_variants.sorted.vcf"
    "out/chr_compress/chr_{chr}_variants.sorted.vcf.gz"
  conda:
    "env/processVCF.yaml"
  log:
    "logs/sortVCF/chr_{chr}.log"
  threads: 4 # just be whatever total n threads is
  shell:
    """ 
    cat {input} | vcf-sort -p {threads} -t out/tmp | \
      bgzip -c > {output}
    """
  


rule split_vcf:
  input:
    "out/all_variants.vcf"
  output:
    "out/chromosomes/chr_{chr}_variants.recode.vcf"
  params:
    prefix = "out/chromosomes/chr_{chr}_variants"
  conda:
    "env/processVCF.yaml"
  log:
    "logs/processVCF/chr_{chr}.log"
  threads: 1
  shell:
    """
    vcftools --vcf {input} --chr {wildcards.chr} --recode --recode-INFO-all --out {params.prefix}
    """


rule call_variants:
  input:
    lst = bam_list
  output:
    "out/all_variants.vcf"
  params:
    lst = lambda wildcards : ",".join(bam_list)
  conda:
    "env/variantPlatypus.yaml"
  log:
    "logs/call_variants.log"
  threads: 99 # just be whatever total n threads is
  shell: 
    """
    platypus callVariants --bamFiles {params.lst} --refFile {config[refFasta]} \
    --filterDuplicates 0 --minMapQual 0 --minFlank 0 --maxReadLength 500 \
    --minGoodQualBases 10 --minBaseQual 20 --nCPU={threads} -o {output}
    """


rule prep_reads:
  input:
    os.path.join(config['bamDir'], "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
  output:
    "out/bam/{sample}_Aligned.sortedByCoord.out.Opossum.bam"
  conda:
    "env/variantPlatypus.yaml"
  log:
    "logs/prep_reads/{sample}.log"
  threads: 4
  shell:
    """ 
    JOBS=$(expr {threads} / 2)
    samtools calmd -b --threads $(($JOBS>1?$JOBS:1)) {input} {config[refFasta]} | \
    python {config[OpossumDir]}/Opossum.py --BamFile /dev/stdin --OutFile {output}
    """


# This rule is not properly implemented / tested yet :~)
# it WILL fail so make sure that the sorted bam files already exist!
# rule align_reads:
#   input: 
#     read1 = ,
#     read2 = ,
#     starGenome, 
#   output:
#     os.path.join("out", "bam", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
#   params:
#     prefix="star/{sample}/{sample}_"
#   conda:
#     "env/alignStar.yaml"
#   log:
#     logs/align_reads/{sample}.log"
#   threads: 99
#   shell:
#     """
#     STAR --twopassMode Basic --genomeDir {input.starGenome} \
#         --readFilesIn {input.read1} {input.read2} \
#         --readFilesCommand zcat --outFileNamePrefix {params.prefix}  \
#         --outSAMtype BAM Unsorted SortedByCoordinate --runThreadN {threads}
#     """


