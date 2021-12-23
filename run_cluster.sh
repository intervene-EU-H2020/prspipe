#!/bin/bash

#SBATCH --job-name=controljob_%j
#SBATCH --output=snakemake_%j.log
#SBATCH --partition=vcpu
#SBATCH --time=48:00:00
#SBATCH -c 1
#SBATCH --mem 4000

# Initialize conda:
eval "$(conda shell.bash hook)"

SNAKEMAKE_ENV='snakemake'

conda activate ${SNAKEMAKE_ENV}

snakemake --snakefile workflow/Snakefile \
          --configfile config/config.yaml \
          --profile ./slurm \
          --directory "${PWD}" \
          "${@}"

