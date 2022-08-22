#!/bin/bash

# script that predicts a given score file

#SBATCH --job-name=sscorer_%j
#SBATCH --output=sscorer_%j.log
#SBATCH --partition=hpcpu
#SBATCH --time=6:00:00
#SBATCH -c 1
#SBATCH --mem 32000

#SBATCH --container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh
#SBATCH --no-container-mount-home
#SBATCH --container-mounts=/dhc/:/dhc/,/home/scratch/:/home/scratch/,/etc/slurm/:/etc/slurm/

if [ $# -ne 4 ]; then
    echo "Need 4 arguments: (1) score-file, (2) biobank, (3) ancestry, (4) output prefix"
    exit 1
fi

Rscript workflow/scripts/R/simple_per_chromosome_scorer.R --target_chr "custom_input/${2}/genotypes/chr" --score $1 --out $4 --keep_pattern '.*_All' --keep "results/${2}/Ancestry_identifier/outlier_detection/AllAncestry.QC.${3}.keep" --afreq "resources/1kg/freq_files/${3}/1KGPhase3.w_hm3.${3}.chr" --mem 32000
