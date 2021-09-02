#!/bin/bash

# RM: TODO - which of these things could be moved to rules instead?

Impute2_1KG_dir="resources/Impute2_1KG"
HapMap3_snplist_dir="resources/HapMap3_snplist"
Geno_1KG_dir="resources/Geno_1KG"

# Create directory for the data
mkdir -p ${Impute2_1KG_dir}

set -e

# Download data using wget
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


# Create directory for the data
mkdir -p ${HapMap3_snplist_dir}
(
cd ${HapMap3_snplist_dir}
if [ ! -f w_hm3.snplist ]; then
    wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
    bunzip2 w_hm3.snplist.bz2
fi
)


# Create directory for the data
mkdir -p ${Geno_1KG_dir}
(
# Download the populations file
cd ${Geno_1KG_dir}
if [ ! -f integrated_call_samples_v3.20130502.ALL.panel ]; then
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
    # Create a version of population file that only contains sample ID, pop and super_pop
fi
cut -f 1-3 integrated_call_samples_v3.20130502.ALL.panel > integrated_call_samples_v3.20130502.ALL.panel_small
#TODO: Create a keep file listing each population super population from the reference.
)

