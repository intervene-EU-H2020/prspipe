#!/bin/bash

#SBATCH --job-name=controljob_%j
#SBATCH --output=snakemake_%j.log
#SBATCH --partition=longrun
#SBATCH --time=48:00:00
#SBATCH -c 1
#SBATCH --mem 3500

SNAKEMAKE_ENV='snakemake'

# Initialize conda:
eval "$(conda shell.bash hook)"
# activate snakemake environment 
conda activate ${SNAKEMAKE_ENV}

snakemake --snakefile workflow/Snakefile \
          --profile ./slurm \
          --directory "${PWD}" \
          "${@}"

