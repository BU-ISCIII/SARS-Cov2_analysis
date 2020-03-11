<img src="./BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

# Artic pipeline for Nanopore data from SARS-Cov2
## RAMPART
### Introduction
First of all we will install RAMPART program. This program will run at the same time of the MinKnow software. It will be displaying the Real-Time coverage of the base-calling, so we can know when to stop the MinIon after having a reliable coverage value.
***_Credit_: https://artic-network.github.io/rampart/***

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
  ***This will take some time as it is going to install all the programs and dependencies required***
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
  ***In our Nanopore machine is installed in /opt/miniconda3/envs/artic-ncov2019 so every user can run it through conda activate /opt/miniconda3/envs/artic-ncov2019***
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
      ***Ths file is not in the repo***
    * `run_configuration.json`: contains information about the current run.
      ***This file is not in the repo***
* Then we have to start a browser introducing in the URL box: _http://localhost:3000/_
