#!/bin/bash

# script that allows exporting the MultiPRS to plink scoring files

#SBATCH --job-name=export_multiprs_%j
#SBATCH --output=export_multiprs_%j.log
#SBATCH --partition=vcpu,hpcpu
#SBATCH --time=6:00:00
#SBATCH -c 1
#SBATCH --mem 32000

#SBATCH --container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh
#SBATCH --no-container-mount-home
#SBATCH --container-mounts=/dhc/:/dhc/,/home/scratch/:/home/scratch/,/etc/slurm/:/etc/slurm/

MODELS=$(find results/ukbb/PRS_evaluation/ -name "*EUR*model.rds")

for m in $MODELS; do
    echo "processing $m"
    Rscript workflow/scripts/R/export_multiPRS.R --model_rds $m
    echo ""
done


