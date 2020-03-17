<img src="../BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

# BU-ISCIII pipeline for  Sequence-Independent, Single-Primer Amplification (SISPA) + Illumina sequencing from SARS-Cov2
## Introduction

## Methodology
The idea is use this pipeline description to create a nextflow pipeline which is in progress, but for the moment we follow a launch mechanism that allows us to structure data and results in a understandable and comprenhensive way.
You can find a description about how we launch or "semi-automatic" pipelines [here](https://github.com/BU-ISCIII/BU-ISCIII/wiki/Execution-guidelines)

## Software and dependencies installation

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
