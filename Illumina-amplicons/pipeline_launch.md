<img src="../BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

# BU-ISCIII pipeline for  Illumina Amplicon sequencing from SARS-Cov2
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
## Pipeline steps:
Amplicon analysis is a little bit different to SISPA analysis because we have to trim the amplicon primers. For the De Novo assembly we are going to trim the primers through sequence. For the consensus steps we are going to trim the primers by position.

### 1. Preprocessing
#### FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc) is used to obtain general quality metrics about the raw reads. It provides information about the quality score distribution across the reads, the per base sequence content (%T/A/G/C), adapter contamination and other overrepresented sequences.

We use our [lablog](./01-fastQC/lablog) as previously explained.
```
bash lablog
```
Running this we obtain the following scripts:
_01_rawfastqc.sh: which performs the quality control of the raw reads:
```
mkdir {sample_id}
fastqc -o {sample_id} --nogroup -t 8 -k 8 ../../00-reads/{sample_id}_R1.fastq.gz ../../00-reads/{sample_id}_R2.fastq.gz
```
_01_unzip.sh: to unzip FastQC results:
```
cd {sample_id}; unzip \*.zip; cd ..
```

#### Trimmomatic
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used to remove adapter contamination and to trim low quality regions. Parameters included for trimming are:
* Nucleotides with phred quality < 10 in 3' end.
* Mean phred quality < 20 in a 4 nucleotide window.
* Read length < 50.

We are going to perform two different trimming steps, one for quality trimming and another one to remove amplicon's primers. The file with the primers sequence is [nCoV-2019.artic.primers.fasta](../data_resources/nCoV-2019.artic.primers.fasta). Thus we are having two different folders.

We move to the quality trimming folder (notrimmedprimers) and we use our [lablog](./02-preprocessing/notrimmedprimers/lablog) as usual:
```
cd notrimmedprimers
bash lablog
```
Then, we obtain the following scripts:
_01_preprocess.sh: To perform the quality trimming of the raw data:
```
mkdir {sample_id}
java -jar /path/to/Trimmomatic/trimmomatic-0.33.jar PE -threads 10 -phred33 ../../00-reads/{sample_id}_R1.fastq.gz ../../00-reads/{sample_id}_R2.fastq.gz {sample_id}/{sample_id}_R1_filtered.fastq {sample_id}/{sample_id}_R1_unpaired.fastq {sample_id}/{sample_id}_R2_filtered.fastq {sample_id}/{sample_id}_R2_unpaired.fastq ILLUMINACLIP:/path/to/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
```
And _02_pgzip.sh: To zip the trimmed fastq files:
```
find . -name "*fastq" -exec pigz -p 5 {} \;
```

Then we move to the amplicon's primer trimming folder (trimmedprimers) and run the [lablog](./02-preprocessing/trimmedprimers/lablog):
```
cd trimmedprimers
bash lablog
```
This creates the following scripts:

_01_preprocess.sh: To perform the primers trimming of the raw data:
```
java -jar /path/to/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 10 -phred33 ../../00-reads/"{sample_id}"_R1.fastq.gz ../../00-reads/"{sample_id}"_R2.fastq.gz {sample_id}/"{sample_id}"_R1_filtered.fastq {sample_id}/"{sample_id}"_R1_unpaired.fastq {sample_id}/"{sample_id}"_R2_filtered.fastq {sample_id}/"{sample_id}"_R2_unpaired.fastq ILLUMINACLIP:../REFERENCES/nCoV-2019.artic.primers.fasta:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
```
And _02_pgzip.sh: To zip the trimmed fastq files:
```
find . -name "*fastq" -exec pigz -p 5 {} \;
```

#### FastQC of the trimmed reads
After the trimmomatic, we perform a fastqc process again to check the quality of the trimmed reads. We are going to perform this quality step for the two trimming steps, thus having two quality folders.

We move to the non trimmed primers folder and use our [lablog](./03-preprocQC/notrimmedprimers/lablog) as previously done:
```
cd notrimmedprimers
bash lablog
```
This created two scripts:
_01_trimfastqc.sh: which performs the quality control of the raw reads:
```
mkdir {sample_id}
fastqc -o {sample_id} --nogroup -t 8 -k 8 ../../02-preprocessing/notrimmedprimers/{sample_id}_R1_filtered.fastq.gz ../../02-preprocessing/notrimmedprimers/{sample_id}_R2_filtered.fastq.gz
```
_01_unzip.sh: to unzip FastQC results:
```
cd {sample_id}; unzip \*.zip; cd ..
```

Then we move to the folder for the quality of the primer trimmed reads and run the [lablog](./03-preprocQC/trimmedprimers/lablog):
```
cd trimmedprimers
bash lablog
```
We obtain the following scripts:
_01_trimfastqc.sh: which performs the quality control of the raw reads:
```
mkdir {sample_id}
fastqc -o {sample_id} --nogroup -t 8 -k 8 ../../02-preprocessing/trimmedprimers/{sample_id}_R1_filtered.fastq.gz ../../02-preprocessing/trimmedprimers/{sample_id}_R2_filtered.fastq.gz
```
_01_unzip.sh: to unzip FastQC results:
```
cd {sample_id}; unzip \*.zip; cd ..
```

### 2. Mapping against host
After performing the preliminary quality controls and trimming, we map the trimmed reads against the host's reference genome. In this case we are going to use the human genome hg38 from the UCSC. For the mapping we use [bwa](http://bio-bwa.sourceforge.net/bwa.shtml) or Burrows-Wheeler Aligner, which is designed for mapping low-divergent sequence reads against reference genomes. The result alignment files are further processed with [SAMtools](http://www.htslib.org/doc/samtools.html), from which sam format is converted to bam, sorted and an index .bai is generated. Finally, Samtools flagstats and [PicardStats](https://broadinstitute.github.io/picard/) are used to obtain statistics over the mapping process.

As for the previous steps we are going to perform two mapping steps, one for the non trimmed primers, and another one for the trimmed primers.

We move to the folder of the non trimmed primers and we run the [lablog](./04-mapping_host/notrimmedprimers/lablog)
```
cd notrimmedprimers
bash lablog
```
From which the following scripts are generated:
_00_mapping.sh: Which is going to perform the mapping of the trimmed reads against the host reference genome:
```
mkdir {sample_id}
bwa mem -t 10 /path/to/host/reference/genome/hg38.fullAnalysisSet.fa ../02-preprocessing/{sample_id}/notrimmedprimers/{sample_id}_R1_filtered.fastq.gz ../02-preprocessing/notrimmedprimers/{sample_id}/{sample_id}_R2_filtered.fastq.gz > {sample_id}/{sample_id}.sam
samtools view -b {sample_id}/{sample_id}.sam > {sample_id}/{sample_id}.bam
samtools sort -o {sample_id}/{sample_id}_sorted.bam -O bam -T {sample_id}/{sample_id} {sample_id}/{sample_id}.bam
samtools index {sample_id}/{sample_id}_sorted.bam
```
_01_flagstat.sh: which is going to perform stats of the mapping through samtools.
```
samtools flagstat {sample_id}/{sample_id}_sorted.bam
```
_02_picadStats.sh: Is going to perform stats about the mapping through Picard.
```
java -jar /path/to/picard-tools-1.140/picard.jar CollectWgsMetrics COVERAGE_CAP=1000000 I={sample_id}/{sample_id}_sorted.bam O={sample_id}/{sample_id}.stats R=/path/to/host/reference/genome/hg38.fullAnalysisSet.fa
```

Then we move to the folder of the trimmed primers and we run the [lablog](./04-mapping_host/trimmedprimers/lablog)
```
cd trimmedprimers
bash lablog
```
From which the following scripts are generated:
_00_mapping.sh: Which is going to perform the mapping of the trimmed reads against the host reference genome:
```
mkdir {sample_id}
bwa mem -t 10 /path/to/host/reference/genome/hg38.fullAnalysisSet.fa ../02-preprocessing/{sample_id}/trimmedprimers/{sample_id}_R1_filtered.fastq.gz ../02-preprocessing/trimmedprimers/{sample_id}/{sample_id}_R2_filtered.fastq.gz > {sample_id}/{sample_id}.sam
samtools view -b {sample_id}/{sample_id}.sam > {sample_id}/{sample_id}.bam
samtools sort -o {sample_id}/{sample_id}_sorted.bam -O bam -T {sample_id}/{sample_id} {sample_id}/{sample_id}.bam
samtools index {sample_id}/{sample_id}_sorted.bam
```
_01_flagstat.sh: which is going to perform stats of the mapping through samtools.
```
samtools flagstat {sample_id}/{sample_id}_sorted.bam
```
_02_picadStats.sh: Is going to perform stats about the mapping through Picard.
```
java -jar /path/to/picard-tools-1.140/picard.jar CollectWgsMetrics COVERAGE_CAP=1000000 I={sample_id}/{sample_id}_sorted.bam O={sample_id}/{sample_id}.stats R=/path/to/host/reference/genome/hg38.fullAnalysisSet.fa
```


### 3.1. Mapping against virus
Once we have mapped the samples to the host, we are going to map the trimmed reads to the reference viral genome. In this care we are going to use the NC_045512.2, Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome. We are going to use the same programs we used for the mapping to the host.

In this case we are just interested in map the non primer trimmed reads to the virus.

We run the [lablog](./05-mapping_virus/lablog)
```
bash lablog
```
We will obtain the following scripts:
_00_mapping.sh: Which is going to perform the mapping of the trimmed reads against the viral reference genome:
```
mkdir {sample_id}
bwa mem -t 10 ../../../REFERENCES/NC_045512.2.fasta ../02-preprocessing/notrimprimers/{sample_id}/{sample_id}_R1_filtered.fastq.gz ../02-preprocessing/notrimprimers/{sample_id}/{sample_id}_R2_filtered.fastq.gz" > {sample_id}/{sample_id}.sam
samtools view -b {sample_id}/{sample_id}.sam > {sample_id}/{sample_id}.bam
samtools sort -o {sample_id}/{sample_id}_sorted.bam -O bam -T {sample_id}/{sample_id} {sample_id}/{sample_id}.bam
samtools index {sample_id}/{sample_id}_sorted.bam
```
_01_flagstat.sh: which is going to perform stats of the mapping through samtools.
```
samtools flagstat {sample_id}/{sample_id}_sorted.bam
```
_02_picadStats.sh: Is going to perform stats about the mapping through Picard.
```
java -jar /path/to/picard-tools-1.140/picard.jar CollectWgsMetrics COVERAGE_CAP=1000000 I={sample_id}/{sample_id}_sorted.bam O={sample_id}/{sample_id}.stats R=../../../REFERENCES/NC_045512.2.fasta
```

### 3.2. Primers positional trimming with iVar.
[iVar](http://gensoft.pasteur.fr/docs/ivar/1.0/manualpage.html) uses primer positions supplied in a BED file to soft clip primer sequences from an aligned, sorted and indexed BAM file. iVar is used to positionally remove the primers of the non trimmed primer sequences.

We run the [lablog](./051-trimPrimers/lablog)
```
bash lablog
```

This creates the follwing scripts:
 _00_ivartrim.sh: Primers trimming
 ```
mkdir -p {sample_id}
samtools view -b -F 4 ../06-mapping_virus/notrimmedprimer/{sample_id}/{sample_id}_sorted.bam > {sample_id}/{sample_id}_onlymapped.bam
samtools index {sample_id}/{sample_id}_onlymapped.bam
ivar trim -e -i {sample_id}/{sample_id}_onlymapped.bam -b ../../../REFERENCES/nCoV-2019.schemeMod.bed -p {sample_id}/{sample_id}_primertrimmed -q 15 -m 50 -s 4
rm {sample_id}/{sample_id}_onlymapped.bam
samtools sort -o {sample_id}/{sample_id}_primertrimmed_sorted.bam -O bam -T {sample_id}/{sample_id} {sample_id}/{sample_id}"_primertrimmed.bam
samtools index {sample_id}/{sample_id}_primertrimmed_sorted.bam
 ```
 _01_flagstat.sh: which is going to perform stats of the mapping through samtools.
 ```
 samtools flagstat {sample_id}/{sample_id}_primertrimmed_sorted.bam
 ```
 _02_picadStats.sh: Is going to perform stats about the mapping through Picard.
 ```
 java -jar /path/to/picard-tools-1.140/picard.jar CollectWgsMetrics COVERAGE_CAP=1000000 I={sample_id}/{sample_id}_primertrimmed_sorted.bam O={sample_id}/{sample_id}.stats R=../../../REFERENCES/NC_045512.2.fasta
 ```


### 4. Variant calling: low freq and mayority calling.
[VarScan](http://varscan.sourceforge.net/) is used for variant calling and two different calls are made:
1. Low frequency variants: we ran VarScan allowing until 3% of alternate allele frequency in order to search for intrahost viral population. This data is going to be meaninful if we have depth of coverages > 1000.
2. Mayority allele calling: we need to generate a consensus genome, so we are going to call just the variants which are present in the mayority of the virus in the sample, we are going to use > 80% in order to keep also some indels that are more difficult called.

We use our [lablog](./06-variant_calling/lablog) as usual, you'd have to change the header for your HPC queue system and parameter:
```Bash
bash lablog
```
Now we are going to get three scripts:
_00_mpileup.sh: varscan requires the generation of a mpileup file. This script will run this command for all the samples:
```
mkdir -p {sample_id}
samtools mpileup -A -d 20000 -Q 0 -f ../../../REFERENCES/NC_045512.2.fasta ../051-trimPrimers/{sample_id}/{sample_id}_primertrimmed_sorted.bam > {sample_id}/{sample_id}.pileup
```
_01_varscan.sh: varscan for low freq variants calling.
```
varscan mpileup2cns {sample_id}/{sample_id}.pileup --min-var-freq 0.02 --p-value 0.99 --variants --output-vcf 1 > {sample_id}/{sample_id}.vcf
```
_02_varscanMajority.sh
```
varscan mpileup2cns {sample_id}/{sample_id}.pileup --min-var-freq 0.8 --p-value 0.05 --variants --output-vcf 1 > {sample_id}/{sample_id}_majority.vcf
```

### 5. Variant effect annotation
We are going to use [SnpEff](http://snpeff.sourceforge.net/) to annotate the variants called in the previous step. The first step is to create the SARS-Cov2 SnpEff database (included in the lablog but we have to do it by hand):
```
cd /processing_Data/bioinformatics/pipelines/miniconda3/envs/virus_illumina_sispa/share/snpeff-4.3.1t-3
vim SnpEff.config ## Add new genome entry.
  sars-cov-2.genome : SARScov2
mkdir -p data/genomes
mkdir -p data/sars-cov-2
cp /path/to/reference/fasta/GCF_009858895.2_ASM985889v3_genomic.fna /path/to/miniconda3/envs/virus_illumina_sispa/share/snpeff-4.3.1t-3/data/genomes/sars-cov-2.fa
cp /path/to/reference/fasta/GCF_009858895.2_ASM985889v3_genomic.gff /path/to/miniconda3/envs/virus_illumina_sispa/share/snpeff-4.3.1t-3/data/data/sars-cov-2/genes.gff
java -jar snpEff.jar build -gff3 -v sars-cov-2
cd /path/to/service/07-annotation/
```
We are going to annotate both low frequency variants and majority variants.
Now we can do the same as always and run the [lablog](./07-annotation/lablog):
```
bash lablog
```
We are going to obtain the following script:
_00_snpEff.sh: To annotate the vcf files:
```
mkdir -p {sample_id}
snpEff sars-cov-2 ../06-variant_calling/{sample_id}_majority.vcf > ${sample_id}_majority.ann.vcf
snpEff sars-cov-2 ../06-variant_calling/{sample_id}.vcf > ${sample_id}_lowfreq.ann.vcf

```
_01_snpsift.sh: Create summary tables:
```
SnpSift extractFields -s "," -e "." {sample_id}/{sample_id}_majority.ann.vcf CHROM POS REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" > {sample_id}/{sample_id}_lowfreq.ann.table.txt
SnpSift extractFields -s "," -e "." {sample_id}/{sample_id}_lowfreq.ann.vcf CHROM POS REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" > {sample_id}/{sample_id}_majority.ann.table.txt
```

### 6. Genome sequence consensus
To obtain the genome sequence consensus we are going to merge the called variants we obtained in the step 4 (Variant calling: low freq and mayority calling.) into the viral reference genome to obtain a consensus between our samples and the reference. For this purposewe are going to use [bgzip](http://www.htslib.org/doc/bgzip.html) and [bcftools](http://www.htslib.org/doc/bcftools.html).

We run the [lablog](./08-mapping_consensus/lablog)
```
bash lablog
```

We obtain the following scripts:
_00_bgzipvcf.sh: To zip the vcf file previously obtained:
```
mkdir {sample_id}_NC_045512.2
bgzip -c ../06-variant_calling/{sample_id}/{sample_id}_majority.vcf > {sample_id}_NC_045512.2/{sample_id}_NC_045512.2.vcf.gz
```
_01_bcftools_index.sh: For indexing the ziped vcf file.
```
bcftools index {sample_id}_NC_045512.2/{sample_id}_NC_045512.2.vcf.gz
```
_02_bcftools_consensus.sh: Create the consensus genome
```
cat ../../../REFERENCES/NC_045512.2.fasta | bcftools consensus {sample_id}_NC_045512.2/{sample_id}_NC_045512.2.vcf.gz > {sample_id}_NC_045512.2/{sample_id}_NC_045512.2_consensus.fasta
```
_03_maskLowCov1.sh: Compute genome coverage
```
bedtools genomecov -bga -ibam ../051-trimPrimers/{sample_id}/{sample_id}_primertrimmed_sorted.bam -g ../../../REFERENCES/NC_045512.2.fasta | awk '\$4 < 20' | bedtools merge > ./{sample_id}_NC_045512/{sample_id}_NC_045512_bed4mask.bed
```
_04_maskLowCov2.sh: Mask low coverage regions
```
bedtools maskfasta -fi ./{sample_id}_NC_045512/{sample_id}_NC_045512_consensus.fasta -bed ./{sample_id}_NC_045512/{sample_id}_NC_045512_bed4mask.bed -fo {sample_id}_NC_045512/{sample_id}_NC_045512_consensus_masked.fasta
```

### 7. De novo assembly
We are going to use the reads that did NOT map to the host reference genome to create a De Novo assembly. First we are going to use Samtools and [Bedtools](https://bedtools.readthedocs.io/en/latest/) to transform the bam file with the unmapped reads into fastq files for the assembly. As the most accurate way to do this depends a lot on the viral genome, the number of reads etc., we are going to perform the De Novo assembly in three different ways:
* Usind [Spades](https://kbase.us/applist/apps/kb_SPAdes/run_SPAdes/release?gclid=Cj0KCQjwjcfzBRCHARIsAO-1_OqnS_JLZtmi_CUOht4Vrl12Rc9bjuTFZKopaxNVWAE6RNgkLjuCfLAaAhjeEALw_wcB):
  * Normal mode
  * Meta Spades mode
* Using [Unicycler](https://github.com/rrwick/Unicycler)
Finally, we are going to use [Quast](http://bioinf.spbau.ru/quast) to determine which assembly is the best for our analysis.

First of all we are going to run the [lablog](./09-assembly/lablog) as allways:
```
bash lablog
```
This will create the following scripts:
_01_unmapped.sh: To select the reads that did not map to the host reference genome.
```
mkdir {sample_id}
samtools view -b -f 4 ../04-mapping_host/{sample_id}/{sample_id}_sorted.bam > {sample_id}/{sample_id}_unmapped.bam
```
_02_name_sorted.sh: To sort the unmapped bam by name:
```
samtools sort -n {sample_id}/{sample_id}_unmapped.bam -o {sample_id}/{sample_id}_unmapped_qsorted.bam
```
_03_bamtofastq.sh: To transofrm the unmapped sorted bam file to fastq files. If the reads are pair-end you can create R1 and R2, if they are single-end you can just create R1.
```
bedtools bamtofastq -i {sample_id}/{sample_id}_unmapped_qsorted.bam -fq {sample_id}/{sample_id}_R1_unmapped.fastq -fq2 {sample_id}/{sample_id}_R2_unmapped.fastq
```
_04_spades_assembly.sh: Normal mode spades de novo assembly.
```
spades.py -t 10 -1 {sample_id}/{sample_id}_R1_unmapped.fastq -2 {sample_id}/{sample_id}_R2_unmapped.fastq -o {sample_id}/{sample_id}
```
_05_spades_assembly_meta.sh: Meta mode spades de novo assembly. This just works for pair-end data. If you have single-end data skip this step.
```
/path/to/spades-3.7.1/bin/spades.py -t 10 -1 {sample_id}/{sample_id}_R1_unmapped.fastq -2 {sample_id}/{sample_id}_R2_unmapped.fastq -o {sample_id}/{sample_id}_meta
```
_06_unicycler.sh: To run unicycler de novo assembly.
```
unicycler -t 10 -o {sample_id} -1 ../02-preprocessing/{sample_id}/{sample_id}_R1_filtered.fastq.gz -2 ../02-preprocessing/{sample_id}/{sample_id}_R2_filtered.fastq.gz
```
Now we are going to rename the output files from the programs:
```
mv {sample_id}/{sample_id}/scaffolds.fasta {sample_id}/{sample_id}/{sample_id}_scaffolds.fasta
mv {sample_id}/{sample_id}_meta/scaffolds.fasta {sample_id}/{sample_id}_meta/{sample_id}_meta_scaffolds.fasta
mv {sample_id}/assembly.fasta {sample_id}/{sample_id}_assembly.fasta
```
_06_quast.sh: To run quast for comparing the different Spades scaffolds:
```
quast.py --output-dir quast_report -R ../../../REFERENCES/NC_045512.2.fasta -G ../../../REFERENCES/NC_045512.2.gff -t 10 $(find . -name "*scaffolds*.fasta" | tr '\n' ' ')"
```
_07_quast.sh: To run quast to compare the unicycler assemblies.
```
quast.py --output-dir quast_report_2 -R ../../../REFERENCES/NC_045512.2.fasta -G ../../../REFERENCES/NC_045512.2.gff -t 10 $(find . -name "*assembly.fasta" | tr '\n' ' ')"
```

### 8. Contig ordering and draft generation.
For reordering the contigs from the spades assembly and generate a draft genome assembly we are going to use [ABACAS](http://abacas.sourceforge.net/index.html).

First we run the [lablog](10-abacas/lablog)
```
bash lablog
```
From which we are going to obtain the following script:
_00_abacas.sh
```
mkdir -p {sample_id}
cd {sample_id}
abacas.pl -r ../../../../REFERENCES/NC_045512.2.fasta -q ../../09-assembly/{sample_id}/{sample_id}/{sample_id}_scaffolds.fasta -m -p nucmer -o {sample_id}_abacas
cd ..
```

### 9. Alignment
We are going to perorm a [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) to align the de novo assemblies to the reference viral genome.

First we run the [lablog](./11-blast/lablog)
```
bash lablog
```
This will create the following scripts:
_00_blast.sh: Performs the alignment.
```
blastn -num_threads 10 -db ../../../REFERENCES/NC_045512.2.fasta -query ../05-assembly/{sample_id}/{sample_id}/{sample_id}_scaffolds.fasta -outfmt \'6 stitle std slen qlen qcovs\' -out {sample_id}_blast.txt
```
_01_filterBlast.sh: Filters the BLAST output to only keep the significant results and parses the results to add headers.
```
awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"}{print \$0,\$5/\$15,\$5/\$14}' {sample_id}_blast.txt | awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"} \$15 > 200 && \$17 > 0.7 && \$1 !~ /phage/ {print \$0}' > {sample_id}_blast_filt.txt; cat header {sample_id}_blast_filt.txt > {sample_id}_blast_filt_header.txt
```
_02_rmTmp.sh: Removes the temporary files created in the previous step.
```
rm *_blast.txt *_blast_filt.txt
```

### 10. Stats and graphs
#### PlasmidID
To graphically see how good or bad the assembly is compared to the reference, we are going to use a program we developped called [plasmidID](https://github.com/BU-ISCIII/plasmidID). Even if the program was created to map assemblies and plasmids, we can use this program to map and align the assemblies to the reference genome and obtain a plot.

We are going to run the [lablog](./12-plasmidID/lablog) as allways:
```
bash lablog
```
From which we are going to obtain the following scripts:
_01_plasmidID_spades.sh: To perform the plasmidID over the normal spades assemblies:
```
plasmidID.sh -d ${directory}/../../../REFERENCES/NC_045512.2.fasta -s {sample_id} -c ${directory}/../05-assembly/{sample_id}/{sample_id}/{sample_id}_scaffolds.fasta -g SPADES --only-reconstruct -C 47 -S 47 -i 60 --no-trim -o ${directory}
```
_02_plasmidID_metaspades.sh: To perorm the plasmidID over the meta spades assemblies:
```
plasmidID.sh -d ${directory}/../../../REFERENCES/NC_045512.2.fasta -s {sample_id} -c ${directory}/../05-assembly/{sample_id}/{sample_id}_meta/{sample_id}_scaffolds_meta.fasta -g META_SPADES --only-reconstruct -C 47 -S 47 -i 60 --no-trim -o ${directory}
```
_03_plasmidID_unicycler.sh: To perform plasmidID over the unicycler assemblies:
```
plasmidID.sh -d ${directory}/../../../REFERENCES/NC_045512.2.fasta -s {sample_id} -c ${directory}/../05-assembly/{sample_id}/{sample_id}_assembly.fasta -g UNICYCLER --only-reconstruct -C 47 -S 47 -i 60 --no-trim -o ${directory}
```
