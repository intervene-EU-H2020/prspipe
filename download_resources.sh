#!/bin/bash

Impute2_1KG_dir="resources/Impute2_1KG"
HapMap3_snplist_dir="resources/HapMap3_snplist"
Geno_1KG_dir="resources/Geno_1KG"
LD_ref_dir="resources/LD_ref"
LDBLOCK_VERSION='ac125e47bf7ff3e90be31f278a7b6a61daaba0dc'

# Create directory for the data
mkdir -p ${Impute2_1KG_dir}

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


# Create a directory for the data
mkdir -p ${LD_ref_dir}/EUR
(
# Download the LD ref score for EUR
cd ${LD_ref_dir}/EUR
if [ ! -f 1.l2.ldscore.gz ]; then
    wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
    tar -xf eur_w_ld_chr.tar.bz2 && rm eur_w_ld_chr.tar.bz2
    mv eur_w_ld_chr/* . && rm -r eur_w_ld_chr
fi
)


# Create a directory for the data
mkdir -p ${LD_ref_dir}/EAS
(
# Download the LD ref score for EAS
cd ${LD_ref_dir}/EAS
if [ ! -f 1.l2.ldscore.gz ]; then
    wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eas_ldscores.tar.bz2
    tar -xf eas_ldscores.tar.bz2 && rm eas_ldscores.tar.bz2
    mv eas_ldscores/* . && rm -r eas_ldscores
fi
)


# Download the LD block information for the DBSLMM method
if [ ! -d ./resources/ldetect-data ]; then
    >&2 echo "Downloading LD block information"
    (
    cd resources
    git clone https://bitbucket.org/nygcresearch/ldetect-data.git && cd ldetect-data && git checkout ${LDBLOCK_VERSION}
    )
else
   (
   cd ./resources/ldetect-data && git checkout ${LDBLOCK_VERSION}
   )
fi

