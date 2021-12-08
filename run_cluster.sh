#!/bin/bash

SNAKEMAKE_ENV='snakemake'

#SBATCH --job-name=controljob_%j
#SBATCH --output=snakemake_%j.log
#SBATCH --partition=vcpu
#SBATCH --time=48:00:00
#SBATCH -c 1
#SBATCH --mem 2000

# Initialize conda:
eval "$(conda shell.bash hook)"

if [ ! -d $snakemake_env ]; then
    ./install.sh
fi

conda activate ${SNAKEMAKE_ENV}

snakemake --snakefile Snakefile \
          --configfile conf/config.yaml \
	  --profile ./slurm \
          --directory "${PWD}" \
	  "${@}"


