#!/bin/bash
#SBATCH -J snakeflow
#SBATCH -p unlimitq
#SBATCH --mem=1G
#SBATCH -o snakemake_output_%j.out
#SBATCH -e snakemake_error_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL

source activate snakemake

snakemake --cores 1 --unlock
snakemake --jobs  10 --cluster-config cluster.yaml --cluster "sbatch --mem {cluster.mem} -c {cluster.cpus}" --use-conda
