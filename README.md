# gapsnake

<img src="images/gapsnake.svg" align="right" /> <br/>*gapsnake* is a <u>snakemake</u>-based workflow for the reconstruction of genome-scale metabolic models using <u>gapseq</u>. It is designed for cases where models for several (from a few to thousands) prokaryotic genomes are reconstructed with the same workflow. It aims to combine an easy and quick usage with full customizability of each individual gapseq module.



## Installation

```sh
# Cloning the development version of gapsnake
git clone https://github.com/Waschina/gapsnake
cd gapsnake

# Create and activate a conda environment "gapsnake"
conda env create -n gapsnake --file gapsnake_env.yaml
conda activate gapsnake
ln -sr `pwd`/gapsnake ${CONDA_PREFIX}/bin/
```

##### Setting up cluster execution

Cluster execution of snakemake workflows is facilitated via [snakemake executor plugins](https://snakemake.github.io/snakemake-plugin-catalog/index.html). To install the plugin for SLURM, you can run:

```sh
conda activate gapsnake
pip install snakemake-executor-plugin-slurm
```

##### CPLEX support

This gapsnake workflow installs gapseq with GLPK as the default LP solver by default. If you have IBM’s ILOG CPLEX installed on your system and want to use it as LP solver instead of GLPK, you can follow the following commands:

First, activate the conda environment if it is not already active:
```sh
conda activate gapsnake
```
Next, start R and run the commands:
```R
install.packages("remotes")
remotes::install_github("Waschina/cobrarCPLEX", configure.args="--with-cplex-dir=/path/to/cplex")
```
Make sure to replace `/path/to/cplex` with the path to your IBM ILOG CPLEX installation (e.g., `/opt/IBM/ILOG/CPLEX_Studio2211/cplex`).

## Quick start

##### Local machine

```sh
# The Setup
cd /path/to/project/directory # Usually one dir up from where your genomes are
conda activate gapsnake


# Initialize gapsnake run
gapsnake init genomes/ # or instead of "genomes/" any path the place where your genomes are
gapsnake recon --cores 16 --keep-going # Set to the maximum number of cores you want to use

```

##### Cluster execution 

```sh
# The Setup
cd /path/to/project/directory # Usually one dir up from where your genomes are
conda activate gapsnake


# Initialize gapsnake run
gapsnake init genomes/ # or instead of "genomes/" any path the place where your genomes are
gapsnake recon --executor slurm --keep-going # Set to the maximum number of cores you want to use
```

Note that the example above is for SLURM job scheduling systems. In case you want to (or have to) use a different scheduling system, please read the [Snakemake Plugin Catalog](https://snakemake.github.io/snakemake-plugin-catalog/index.html). 
