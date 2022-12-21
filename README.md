# gapsnake

<img src="images/gapsnake.svg" align="right" /> <br/>*gapsnake* is a <u>snakemake</u>-based workflow for the reconstruction of genome-scale metabolic models using <u>gapseq</u>. It is designed for cases where models for several (from a few to thousands) prokaryotic genomes are reconstructed with the same workflow. It aims to combine an easy and quick usage with full customizability of each individual gapseq module.



## Prerequisites

- you will need to have *gapseq* and its dependencies installed. Also, make sure that you have `gapseq` in your PATH variable.

## Installation

```sh
# Cloning the development version of gapsnake
git clone https://github.com/Waschina/gapsnake
cd gapsnake

# Create and activate a conda environment "gapsnake"
conda env create -n gapsnake --file gapsnake_env.yml
conda activate gapsnake
```

##### Setting up cluster execution

Note: Until now, only tested with the SLURM scheduler on CAU's caucluster.

```sh
conda activate gapsnake
cookiecutter --output-dir ~/.config/snakemake https://github.com/metagenome-atlas/clusterprofile.git
```



## Quick start

Local machine

```R
# The Setup
cd /path/to/project/directory # Usually one dir up from where your genomes are
conda activate gapsnake
export PATH="${PATH}:/path/to/gapsnake/"

# Initialize gapsnake run
gapsnake init genomes/ # or instead of "genomes/" any path the place where your genomes are


```



Cluster execution

```sh
gapsnake recon --profile cluster
```

