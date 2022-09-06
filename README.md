# OBITools workflow

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6676577.svg)](https://doi.org/10.5281/zenodo.6676577)

<!-- TABLE OF CONTENTS -->
<details>
    <summary>Table of Contents</summary>
    <ol>
        <li>
            <a href="#about">About</a>
        </li>
        <li>
            <a href="#getting-started">Getting Started</a>
            <ul>
                <li><a href="#installation">Installation</a></li>
            </ul>
            <ul>
                <li><a href="#directories-and-files-structure">Directories and files structure</a></li>
            </ul>
            <ul>
                <li><a href="#download-your-data">Download your data</a></li>
            </ul>
        </li>
        <li>
            <a href="#usage">Usage</a>
            <ul>
                <li><a href="#configuration">Configuration</a></li>
            </ul>
            <ul>
                <li><a href="#run-the-workflow">Run the workflow</a></li>
        </li>
        <li>
            <a href="#example">Example</a>
            <ul>
                <li><a href="#download-toy-data">Download toy data</a></li>
            </ul>
            <ul>
                <li><a href="#run-the-workflow">Run the workflow</a></li>
            </ul>
            <ul>
                <li><a href="#option-merging-libraries">Option: merging libraries</a></li>
            </ul>
        </li>
        <li>
            <a href="#going-further">Going further</a>
        </li>
        <li>
            <a href="#acknowledgements">Acknowledgements</a>
        </li>
        <li>
            <a href="#how-to-cite-this-repository">How to cite this repository</a>
        </li>
        <li>
            <a href="#references">References</a>
        </li>
    </ol>
</details>

## About

This is a snakemake workflow based on the obitools suite of programs, that analyzes DNA metabarcoding data.

Sequence analysis is performed with the obitools (Boyer et al. 2016) and sumaclust (Mercier et al. 2013) through a Snakemake pipeline (Mölder et al. 2021).


## Getting started

### Installation

#### Dependencies

In order to run the workflow, you must have installed the following programs:

- [python3](https://www.python.org/downloads/)
- [conda](https://docs.conda.io/en/latest/miniconda.html)
- [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

Please note that the workflow is currently running exclusively on Unix systems.

#### Install the workflow

Clone the repository:
```sh
git clone https://github.com/AnneSoBen/obitools_workflow.git
```

### Directories and files structure

The repository contains five folders:
- `config/`: contains the configuration file of the Snakemake workflow (`config.yaml`). This is where the value of the options for the various commands used is defined.
- `log/`: where log files of each rule are written.
- `resources/`: where you should download/copy your raw data (cf. _Download your data_)
- `results/`: where all output files are written.
- `workflow/`: contains the Snakemake workflow (`Snakefile`), the configuration file of the submission parameters on the cluster (`cluster.yaml`), the script to submit the workflow on the cluster (`sub_smk.sh`). 

### Download your data

Download/copy your data in the `resources/` folder. Three files are required:
- forward and reverse fastq files
- the corresponding ngsfilter file

They should be named as follows: `prefix_R1.fastq`, `prefix_R2.fastq`, `prefix_ngsfilter.tab`

And be put in a subfolder whose name is the prefix of the files (see _Example_).


## Usage

### Configuration

Before running the workflow, the configuration file (`config/config.yaml`) has to be edited. The parameters that can be set are listed in the table below:

| parameter          | description                                                                          | concerned rule(s)                                                                                    | default value | comment                                                                                                                                                              |
|--------------------|--------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|---------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| tomerge            | whether to merge libraries before dereplication                                      | merge_demultiplex                                                                                    | FALSE         | should be set to 'TRUE' if you analyse several libraries and that you want to merge them                                                                             |
| resourcesfolder    | relative path to the folder containing resource files (fastq files and ngsfilter)    | split_fastq, demultiplex                                                                             | ../resources  | should not be changed, unless you want to rename the folder                                                                                                          |
| resultsfolder      | relative path to the folder where output files will be written                       | all                                                                                                  | ../results    | should not be changed, unless you want to rename the folder                                                                                                          |
| fastqfiles         | prefix of the name of the resource fastq files and ngsfilter                         | all                                                                                                  | wolf_diet     | must be changed to match your files name prefix                                                                                                                      |
| mergedfile         | prefix of the name of the output files if tomerge=TRUE                               | merge_demultiplex, split_fasta, derepl, merge_derepl, basicfilt, clustering, merge_clust, tab_format | wolf_diet     | must be changed for the merged files name prefix you want                                                                                                            |
| split_fastq:nfiles | number of files to create when splitting fastq files for pairing                     | split_fastq                                                                                          | 2             | should be changed according to the size of you dataset: the bigger it is, the more you will want to split your initial files - useful only on multi-threaded systems |
| minscore           | minimum alignment score required for pairing                                         | alifilt                                                                                              | 40.00         | set according to Taberlet et al. 2018                                                                                                                                |
| split_fasta:nfiles | number of files to create when splitting demultiplexed fasta files for dereplication | split_fasta                                                                                          | 2             | should be changed according to the size of you dataset: the bigger it is, the more you will want to split your initial file(s)                                       |
| minlength          | minimum sequence length (in bp)                                                      | basicfilt                                                                                            | 80            | must be changed according to the minimum length expected for your barcode                                                                                            |
| mincount           | minimum number of reads per unique sequence                                          | basicfilt                                                                                            | 1             | it's up to you!                                                                                                                                                      |
| minsim             | similarity threshold for clustering                                                  | clustering                                                                                           | 0.97          | it's up to you!                                                                                                                                                      |


If you run the workflow on a SLURM cluster, you must also check the `workflow/cluster.yaml` that sets up the ressources available for each rule.

### Run the workflow

Then, run the workflow:
```sh
cd workflow
conda activate snakemake
snakemake -c1 --use-conda
```

Alternatively, you can run the workflow in a single command on a SLURM cluster by submitting the `sub_smk.sh` file:
```sh
cd workflow
sbatch sub_smk.sh
```


## Example

### Download toy data

If you want to test the workflow, download toy data from the obitools tutorial (https://pythonhosted.org/OBITools/wolves.html) in the `resources/` folder:
```sh
wget -O resources/wolf_tutorial.zip https://pythonhosted.org/OBITools/_downloads/wolf_tutorial.zip
unzip resources/wolf_tutorial.zip -d resources/
mv resources/wolf_tutorial resources/wolf_diet
rm resources/wolf_tutorial.zip
```
Rename the files to fit the template decribed above (or create symbolic links):
```sh
ln -s wolf_F.fastq resources/wolf_diet/wolf_diet_R1.fastq
ln -s wolf_R.fastq resources/wolf_diet/wolf_diet_R2.fastq
ln -s wolf_diet_ngsfilter.txt resources/wolf_diet/wolf_diet_ngsfilter.tab
```
You should get this directories and files structure:
```sh
tree
```

```
.
├── config
│   └── config.yaml
├── LICENSE
├── log
├── README.md
├── resources
│   └── wolf_diet
│       ├── db_v05_r117.fasta
│       ├── embl_r117.ndx
│       ├── embl_r117.rdx
│       ├── embl_r117.tdx
│       ├── wolf_diet_ngsfilter.tab -> wolf_diet_ngsfilter.txt
│       ├── wolf_diet_ngsfilter.txt
│       ├── wolf_diet_R1.fastq -> wolf_F.fastq
│       ├── wolf_diet_R2.fastq -> wolf_R.fastq
│       ├── wolf_F.fastq
│       └── wolf_R.fastq
├── results
└── workflow
    ├── cluster.yaml
    ├── Snakefile
    └── sub_smk.sh
```

Note that the name of the subfolder containing your source files (fastq and ngsfilter files) should be the prefix of the files.

The config.yaml file is already modified to fit this data.

### Run the workflow

Now run the workflow:
```sh
cd workflow/
conda activate snakemake
snakemake -c1 --use-conda
```

### Option: merging libraries

You may want to merge libraries, for example if technical replicates are splitted in different libraries. To allow this, the value of "tomerge" in the `config/config.yaml` file should be set to `TRUE`. Besides, the prefix of your library files should be listed in the `config/config.yaml` file, such as:

```
tomerge:
  TRUE
resourcesfolder:
  ../resources/
resultsfolder:
  ../results/
fastqfiles:
  - myfirstlibfileprefix
  - mysecondlibfileprefix
mergedfile:
  mymergedlibs
```

The source files of each library should be in separate subfolders. For example:

```
└─ resources
 └── myfirstlibprefix
 |   ├── myfirstlibprefix_ngsfilter.tab
 |   ├── myfirstlibprefix_R1.fastq
 |   └── myfirstlibprefix_R2.fastq
 └── mysecondlibprefix
     ├── mysecondlibprefix_ngsfilter.tab
     ├── mysecondlibprefix_R1.fastq
     └── mysecondlibprefix_R2.fastq
```

Two ngsfilter files will be necessary: `resources/myfirstlibfileprefix/myfirstlibfileprefix_ngsfilter.tab` and `resources/myfirstlibfileprefix/mysecondlibfileprefix_ngsfilter.tab`.

:warning: If you want to be able to distinguish your technical replicates in the final output, don't forget to give your samples different names in the ngsfilter files, e.g. for a sample named "sample", you could change its name to "sample_a" in the first ngsfilter file and "sample_b" in the second ngsfilter file (if you have to technical replicates).

The value of "mergedfile" corresponds to the prefix of the merged files from the dereplication to the end of the workflow.


## Going further

You may want to clean up potential molecular artifacts: have a look at the R package [metabaR](https://github.com/metabaRfactory/metabaR)!


## Acknowledgements

Thanks to **[Lucie Zinger](https://luciezinger.wordpress.com/)**, **[Frédéric Boyer](https://www.researchgate.net/profile/Frederic-Boyer-3)**, **[Céline Mercier](https://www.celine-mercier.info/)** and **Clément Lionnet** for their help with the obitools! Also thanks to the **[ECOFEED](https://cordis.europa.eu/project/id/817779/fr)** project for funding the development of the first version of this workflow.


## How to cite this repository

Anne-Sophie Benoiston. (2022). AnneSoBen/obitools_workflow: v1.0.2. GitHub. https://doi.org/10.5281/zenodo.6676577.

:triangular_flag_on_post: Don't forget to cite this repository if you use it for your research :slightly_smiling_face:


## References

Boyer, F., Mercier, C., Bonin, A., Bras, Y. L., Taberlet, P., & Coissac, E. (2016). obitools: A unix-inspired software package for DNA metabarcoding. Molecular Ecology Resources, 16(1), 176‑182.

Mercier, C., Boyer, F., Bonin, A., & Coissac, E. (2013). SUMATRA and SUMACLUST: fast and exact comparison and clustering of sequences. In Programs and Abstracts of the SeqBio 2013 workshop. Abstract (pp. 27-29).

Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., ... & Köster, J. (2021). Sustainable data analysis with Snakemake. F1000Research, 10.

Zinger, L., Lionnet, C., Benoiston, A. S., Donald, J., Mercier, C., & Boyer, F. (2021). metabaR: an R package for the evaluation and improvement of DNA metabarcoding data quality. Methods in Ecology and Evolution, 12(4), 586-592.
