<img src="../BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

# BU-ISCIII pipeline for  Sequence-Independent, Single-Primer Amplification (SISPA) + Illumina sequencing from SARS-Cov2
## Introduction

## Methodology
The internal structure of each analysis folder (folder naming: **{date}_{CUSTOMNAME}_{user/institution}**) consist on 6 folders:
1. **ANALYSIS:** each analysis (scripts, documentation and results) we make will be in a folder here, named as XX-name where XX is the number of the analysis and name is the name of the analysis.
2. **DOC:** any documentation required for the service or talks, posters an publications derived from its results.
3. **RAW:** raw input data is copied here while we are working on the service
4. **REFERENCES:** reference sequences and files for the analysis
5. **RESULTS:** one folder per each results submission for the service. They are named as YYYYMMDD_entregaXX, where YYYYMMDD is the date of submition in the numerical format year-month-day, and XX is the number of the submission.This folder will be deprecated soon because this information will be in bioinfo_doc/services
6. **TMP:** temporal files

The idea is use this pipeline description to create a nextflow pipeline which is in progress, but for the moment we follow a mechanism for running the different steps that allows us to structure data and results in a understandable and comprenhensive way.
You can find a description about how we launch or "semi-automatic" pipelines [here](https://github.com/BU-ISCIII/BU-ISCIII/wiki/Execution-guidelines)

## Software and dependencies installation
We use a conda environment which provides all the sofware needed.
## Pipeline steps

### 1. Preprocessing
### 2. Mapping against host
### 3. Mapping against virus
### 4. Variant calling: low freq and mayority calling.
[VarScan]() is used for variant calling and two different calls are made:
1. Low frequency variants: we ran VarScan allowing until 3% of alternate allele frequency in order to search for intrahost viral population. This data is going to be meaninful if we have depth of coverages > 1000.
2. Mayority allele calling: we need to generate a consensus genome, so we are going to call just the variants which are present in the mayority of the virus in the sample, we are going to use > 80% in order to keep also some indels that are more difficult called.

We use our [lablog](./06-variant_calling/lablog) as usual, you'd have to change the header for your HPC queue system and parameter:
```Bash
bash lablog
```
Now we are going to get three scripts:
_00_mpileup.sh: varscan requires the generation of a mpileup file. This script will run this command for all the samples:
```
samtools mpileup -A -d 20000 -Q 0 -f ../../../REFERENCES/NC_045512.2.fasta ../06-mapping_virus/${sample_id}/${sample_id}_sorted.bam > ${sample_id}.pileup
```
_01_varscan.sh: varscan for low freq variants calling.
```
varscan mpileup2cns ./\${sample_id}.pileup --min-var-freq 0.02 --p-value 0.99 --variants --output-vcf 1 > \${sample_id}.vcf
```
_02_varscanMajority.sh
```
varscan mpileup2cns ./\${sample_id}.pileup --min-var-freq 0.8 --p-value 0.05 --variants --output-vcf 1 > \${sample_id}_mayority.vcf
```

### 5. Variant effect annotation
### 5. Genome sequence consensus
### 6. De novo assembly
### 7. Contig ordering and draft generation.
### 8. Stats and graphs
