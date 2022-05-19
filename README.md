# obitools_workflow

## About

This is a snakemake workflow based on the obitools suite of programs, that analyzes DNA metabarcoding data.

Sequence analysis is performed with the obitools (Boyer et al. 2016) and sumaclust (Mercier et al. 2013) through a Snakemake pipeline (Molder et al. 2021).

## Getting started

### Prerequisites

This workflow is meant to be executed on a computing cluster running with **SLURM**. It has been written to run on the Genotoul computing cluster (http://bioinfo.genotoul.fr/).

### Installation

Clone the repository:
```sh
git clone ...
```
### Directories and files structure

The repository contains five folders:
- config: contains the configuration file of the Snakemake workflow (`config.yaml`). This is where the value of the options for the various commands used is defined.
- log: where log files of each rule are written.
- resources: where you should download/copy your raw data (cf. _Download your data_)
- results: where all output files are written.
- workflow: contains the Snakemake workflow (`Snakefile`), the configuration file of the submission parameters on the cluster (`cluster.yaml`), the script to submit the workflow on the cluster (`sub_smk.sh`). 

### Download your data

Download/copy your data in the `resources/` folder. Three files are required:
- forward and reverse fastq files
- the corresponding ngsfilter file

They should be named as follows:
basename_R1.fastq
basename_R2.fastq
basename_ngsfilter.tab

## Usage

Before running the workflow, the two configuration files have to be modified: `workflow/cluster.yaml` that sets up the ressources available for each rule, and `config/config.yaml` where you can edit the values of the parameters used by the rules and the basename of your files.

Then, to run the workflow in a single command on the cluster:

```
sbatch sub_smk.sh
```

## Example

### Download toy data

If you want to test the workflow, download toy data from the obitools tutorial (https://pythonhosted.org/OBITools/wolves.html) in the resources folder:
```sh
wget -O resources/wolf_tutorial.zip https://pythonhosted.org/OBITools/_downloads/wolf_tutorial.zip
unzip resources/wolf_tutorial.zip -d resources/
mv resources/wolf_tutorial/* resources/
rm resources/wolf_tutorial.zip
rm -r resources/wolf_tutorial
```
Rename the files to fit the template decribed above (or create symbolic links):
```sh
ln -s wolf_F.fastq resources/wolf_diet_R1.fastq
ln -s wolf_R.fastq resources/wolf_diet_R2.fastq
ln -s wolf_diet_ngsfilter.txt resources/wolf_diet_ngsfilter.tab
```

The config.yaml file is already modified to fit this data.

## Acknowledgements

Thanks to **Lucie Zinger**, **Frédéric Boyer**, **Céline Mercier** and **Clément Lionnet** for their help with the obitools!
