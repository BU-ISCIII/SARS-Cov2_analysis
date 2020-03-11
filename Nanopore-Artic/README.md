<img src="./BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

# Artic pipeline for Nanopore data from SARS-Cov2
## RAMPART
### Introduction
First of all we will install RAMPART program. This program will run at the same time of the MinKnow software. It will be displaying the Real-Time coverage of the base-calling, so we can know when to stop the MinIon after having a reliable coverage value.

The link for more documentation is [here](https://artic.network/ncov-2019/ncov2019-using-rampart.html)

### Installation
The first thin is setting up the computing environment ad  described in: [ ARTIC-nCoV-ITSetup](https://artic.network/ncov-2019/ncov2019-it-setup.html)

* First we will need to install the ARTIC nCoV-2019 data and software repository:
```
git clone --recursive https://github.com/artic-network/artic-ncov2019.git
```
* Then we will create the conda environment with the yml in the previous repo:
```
conda env create -f artic-ncov2019/environment.yml
```
  ***This will thake come time as it is going to install all the programs and dependencies requiered***
* Then we will test if the installation went right:
***In out Nanopore machine is installed in /opt/miniconda3/envs/artic-ncov2019 so every user can run it through conda activate /opt/miniconda3/envs/artic-ncov2019***
```
conda activate artic-ncov2019
rampart --help
```
