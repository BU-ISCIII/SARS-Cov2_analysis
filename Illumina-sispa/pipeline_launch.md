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

Here you can see an example of the tree of the working folder:

```
20200306_SARS-Cov2_RESPIRATORIOS_IC_C
├── ANALYSIS
│   ├── 00-reads
│   ├── 20200313_illumina_sispa
│   └── 20200313_nanopore
├── DOC
│   ├── galaxy_workflow_cov19
│   └── GalaxyWorkflow_SARS-CoV-2.sh
├── RAW
│   ├── logs
│   ├── QC
│   ├── SRR11140744_1.fastq.gz
│   ├── SRR11140744_2.fastq.gz
│   ├── SRR11140746_1.fastq.gz
│   ├── SRR11140746_2.fastq.gz
│   ├── SRR11140748_1.fastq.gz
│   ├── SRR11140748_2.fastq.gz
│   ├── SRR11140750_1.fastq.gz
│   └── SRR11140750_2.fastq.gz
├── REFERENCES
│   ├── GCF_009858895.2_ASM985889v3_genomic.fna
│   ├── GCF_009858895.2_ASM985889v3_genomic.gff
│   ├── NC_045512.2.fasta -> GCF_009858895.2_ASM985889v3_genomic.fna
│   ├── NC_045512.2.fasta.amb
│   ├── NC_045512.2.fasta.ann
│   ├── NC_045512.2.fasta.blast.tmp.nhr
│   ├── NC_045512.2.fasta.blast.tmp.nin
│   ├── NC_045512.2.fasta.blast.tmp.nsq
│   ├── NC_045512.2.fasta.bwt
│   ├── NC_045512.2.fasta.fai
│   ├── NC_045512.2.fasta.nhr
│   ├── NC_045512.2.fasta.nin
│   ├── NC_045512.2.fasta.nog
│   ├── NC_045512.2.fasta.nsd
│   ├── NC_045512.2.fasta.nsi
│   ├── NC_045512.2.fasta.nsq
│   ├── NC_045512.2.fasta.pac
│   ├── NC_045512.2.fasta.sa
│   ├── NC_045512.2.gff -> GCF_009858895.2_ASM985889v3_genomic.gff
├── RESULTS
└── TMP
    ├── ZEBOV_3Samples_NB
    └── ZEBOV_3Samples_NB_MinIT_guppy.tgz
```

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
