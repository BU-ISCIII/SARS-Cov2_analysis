<img src="../BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

# Artic pipeline for Nanopore data from SARS-Cov2
## RAMPART
### Introduction
First of all we will install RAMPART program. This program will run at the same time of the MinKnow software. It will be displaying the Real-Time coverage of the base-calling, so we can know when to stop the MinIon after having a reliable coverage value.
>_Credit_: https://artic-network.github.io/rampart/

Link for more documentation [RAMPART](https://artic.network/ncov-2019/ncov2019-using-rampart.html)

### Installation
This installation protocol is specific for the 2019 Cov. The first thing is setting up the computing environment as described in: [ ARTIC-nCoV-ITSetup](https://artic.network/ncov-2019/ncov2019-it-setup.html)

* First we will need to install the ARTIC nCoV-2019 data and software repository:
```
git clone --recursive https://github.com/artic-network/artic-ncov2019.git
```
* Then we will create the conda environment with the .yml in the previous repo:
```
conda env create -f artic-ncov2019/environment.yml
```
  >This will take some time as it is going to install all the programs and dependencies required

  * We found an errore when we tried to install the environment in the cluster due to compiler errors. We solved it doing the following:
  ```
  module load Compilers/gcc-5.5.0
  module load Libraries/glibc-2.14
  CXX=/opt/Compilers/gcc-5.5.0/bin/g++
  CC=/opt/Compilers/gcc-5.5.0/bin/gcc
  conda install -c conda-forge binutils
  pip install git+https://github.com/artic-network/Porechop.git@v0.3.2pre
  ```

* Then we will test if the installation went right:
  >In our Nanopore machine is installed in /opt/miniconda3/envs/artic-ncov2019 so every user can run it through conda activate /opt/miniconda3/envs/artic-ncov2019
```
conda activate artic-ncov2019
rampart --help
```

### Running RAMPART
* Prepare the sequencing run with MinKnow by creating a specific directory for the run. The name of the folder should be the same as the Experiment name introduced in MinKnow ('Experiment' tab in the left panel).
```
mkdir 01-RAMPARD/<dir_name>
cd 01-RAMPARD/<dir_name>
```
* For RAMPARD to work properly we should set MinKnow setting to write just 1000 reads per FASTQ file ('Output' tab in the left panel).
  * File type = FAST5 / Reads per file =4000.
  * File type = FASTQ / Reads per file = 1000.
  * Remember the output location path of your files <out_dir>.
* If you are using ONT native barcoding to multiplex samples you can optionally create _barcodes.csv_. This file should have one line per barcode and only contain the barcodes that are really contained in the library.
* Whenever the fist reads are basecalled we can run RAMPART:
```
rampart --protocol <path_to_repositories>/artic-ncov2019/rampart --basecalledPath <out_dir>/<dir_name>/fastq_pass
```
  * `--protocol` defined the path where the configuration files are stored. Normally the protocol directory is virus-specific, not run-specific. The one we have in the repository is specific for the novel 2019 Cov. [We can modify the configuration files for our own run](https://artic-network.github.io/rampart/docs/protocols):
    * `protocol.json`: Description of the protocol’s purpose
    * `genome.json`: describes the reference genome of what’s being sequenced.
    * `primers.json`: describes the position of primers across the genome. This file will be used by RAMPART to draw the the amplicons in the coverage plots. If it is not present then no amplicons will be shown in RAMPART.
    * `pipelines.json`: describes the pipelines used for data processing and analysis by RAMPART.
      >Ths file is not in the repo

    * `run_configuration.json`: contains information about the current run.
      >This file is not in the repo

* Then we have to start a browser introducing in the URL box: _http://localhost:3000/_
## MinIon and MinKnow
We will set all the MinKnow options as described in [this protocol](https://www.protocols.io/view/ncov-2019-sequencing-protocol-bbmuik6w).
* Choose flow cell 'FLO-MIN106' from the drop-down menu.
* Then select the flowcell so a tick appears.
* Click the 'New Experiment' button in the bottom left of the screen.
* On the New experiment popup screen, select the running parameters for your experiment from the individual tabs:
  * Experiment:
    * Name the run in the experiment field the same as the directory we created for RAMPART.
    * Leave the sample field blank.
  * Kit:
    * Selection: Select LSK109 as there is no option for native barcoding (NBD104).
  * Run Options:
    * Set the run length to 6 hours (you can stop the run once sufficient data has been collected as determined using RAMPART).
  * Basecalling:
    * Leave basecalling turned but select 'fast basecalling'.
  * Output:
    * Remember the location directory.
    * Reduce Reads per file for FASTQ to 1000.
  * Start run.

## Artic Bioinformatic pipeline
### Create directory for the analysis
```
mkdir ANALYSIS
cd ANALYSIS
mkdir 02-artic_pipeline
```
### Activate the ARTIC environment
```
conda activate artic-ncov2019
```
> We won't run Basecalling with Guppy as we will do it with MinKnow.

### Consensus sequence generation
This will collect all the FASTQ files that have been created by the basecaller into a single file. We will filter the length of the reads to the ones that are between 400 and 700 to remove chimeric reads.
> We can use the minimum lengths of the amplicons as the minim length and the maximum length of the amplicons plus 200 as the maximum. I.e. if your amplicons are 300 base pairs, use –min-length 300 –max-length 500.
```
artic gather --min-length <amplicon_length> --max-length <amplicon_length+200> --prefix <prefix_gathered_files> --directory <path_to_basecalled_reads>
```
> We recomend to put the prefix_gathered_files the name as the run_name because it will be the output fastq file.

### Demultiplexing with Porechop
This step is mandatory.
```
artic demultiplex --threads 4 <prefix_gathered_files>.fastq
```
> This is running porechop as: porechop --verbosity 2 --untrimmed -i \"%s\" -b %s --native_barcodes --discard_middle --require_two_barcodes --barcode_threshold 80 --threads %s --check_reads 10000 --barcode_diff 5 > %s.demultiplexreport.txt" % (args.fasta, tmpdir, args.threads, args.fasta))
  > -verbosity 2: shows the actual trimmed/split sequences for each read (described more below).
  > --untrimmed: Bin reads but do not trim them
  > -i: imput reads
  > -b: Reads will be binned based on their barcode and saved to separate files in this directory (incompatible with --output). The program creates a temporary directory to store this files and then removes it after renameming the output files.
  > --native_barcodes: Only attempts to match the 24 native barcodes.
  > --discard_middle: Reads with middle adapters will be discarded (default: reads with middle adapters are split) (required for reads to be used with Nanopolish).
  > --require_two_barcodes: Reads will only be put in barcode bins if they have a strong match for the barcode on both their start and end (default: a read can be binned with a match at its start or end).
  > --barcode_threshold: A read must have at least this percent identity to a barcode to be binned (default: 75.0).
  > --check_reads: This many reads will be aligned to all possible adapters to determine which adapter sets are present (default: 10000).
  > --barcode_diff: If the difference between a read's best barcode identity and its second-best barcode identity is less than this value, it will not be put in a barcode bin (to exclude cases which are too close to call)(default: 5.0).

### Create nanopolish index
We only have to do this once per sequencing run.
```
nanopolish index -s run_name_sequencing_summary.txt -d <path_to_fast5> <prefix_gathered_files>.fastq
```
> nanopolish index: Build an index mapping from basecalled reads to the signals measured by the sequencer
> -s: the sequencing summary file from albacore, providing this option will make indexing much faster.
> -d path to the directory containing the raw ONT signal files. This option can be given multiple times.

### Run the MinION pipeline
We have to run this command once for each of the barcode file.
```
artic minion --normalise 200 --threads 4 --scheme-directory ~/artic-ncov2019/primer-schemes --read-file run_name_pass_<barcode_ID>.fastq --nanopolish-read-file <prefix_gathered_files>.fastq nCoV-2019/V1 <samplename_for_barcode_ID>
```
> --normalise: Normalise down to moderate coverage to save runtime.
> --scheme-directory: path to scheme-directory that is included in the github repository. This folder contains:
> --read-file: path to the fastq file demultiplexed by Porechop for ONE barcode.
> --nanopolish-read-file: path to the gathered fastq file thas has been indexed with nanopolish.
> nCoV-2019/V1: path to the version of the primer-schemes we are going to use.
  > nCoV-2019.log
  > nCoV-2019.pdf
  > nCoV-2019.pickle
  > nCoV-2019.reference.fasta: MN908947.3 (Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome)
  > nCoV-2019.scheme.bed
  > nCoV-2019.svg
  > nCoV-2019.tsv
  > nCoV-2019_SMARTplex.tsv
> <samplename_for_barcode_ID>: name of the sample that is identified with the barcode_ID.

**Output files:**
* <samplename_for_barcode_ID>.primertrimmed.bam
  * BAM file for visualisation after primer-binding site trimming
* <samplename_for_barcode_ID>.vcf
  * Detected variants in VCF format.
* <samplename_for_barcode_ID>.variants.tab
  * Detected variants.
* <samplename_for_barcode_ID>.consensus.fasta
  * Consensus sequence


If we would want to put all the consensus sequences in one file:
```
cat *.consensus.fasta > <output_consensus_genomes>.fasta
```

### Visualize genomes in Tablet
* Open a new terminal:
  ```
  conda activate tablet
  tablet
  ```
* Go to "Open Assembly"
  * Load the BAM as the first file
  * Load the reference file as the second file: artic-ncov2019/primer_schemes/nCoV-2019/V1/nCoV-2019.reference.fasta (MN908947.3)
  * Select variants mode in Color Schemes for ease of viewing variants.
