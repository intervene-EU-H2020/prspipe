#!/bin/bash

# RM: TODO - get rid of this script and replace with rules.
# currenlty I think the Impute2 data is not used at all...
# other things that used to be here have been moved to rules in setup.smk

Impute2_1KG_dir="resources/Impute2_1KG"

# Create directory for the data
mkdir -p ${Impute2_1KG_dir}

set -e

(
cd ${Impute2_1KG_dir}
if [ ! -f 1000GP_Phase3.tgz ]; then
    wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz
    tar -zxvf 1000GP_Phase3.tgz
fi
if [ ! -f 1000GP_Phase3_chrX.tgz ]; then
    wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz
    tar -zxvf 1000GP_Phase3_chrX.tgz
fi
)
