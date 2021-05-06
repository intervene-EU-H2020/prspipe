#!/bin/bash

# Date:         24.03.2021
# Author:       Remo Monti & Sophie Wharrie
# Description:  Download binaries, move them to bin, clone GenoPred repo

GENOPRED_VERION="latest"
LDSC_VERSION="aa33296abac9569a6422ee6ba7eb4b902422cc74"
PRSICE_VERSION="2.3.3"

if [ ! -f README.md ]; then
    echo "Error: Are you in the project base directory? Abort."
    exit 1
fi

if [ ! -d ./bin ]; then
    mkdir bin
fi

# Install Plink 2.0s
if [ ! -f ./bin/plink2 ]; then
    >&2 echo "Downloading Plink 2.0 binaries"
    (
    cd bin
    wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_avx2.zip
    unzip plink2_linux_avx2.zip && rm plink2_linux_avx2.zip
    )
fi

# Install Plink 1.9
if [ ! -f ./bin/plink ]; then
    >&2 echo "Downloading Plink 1.9 binaries"
    (
    cd bin
    wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip
    unzip plink_linux_x86_64_20201019.zip && rm plink_linux_x86_64_20201019.zip
    )
fi

# Install GCTA
# TODO: add links to the binaries, or move binaries to ./bin
if [ ! -f ./bin/gcta/gcta64 ]; then
    >&2 echo "Downloading GCTA 1.93.2 binaries"
    (
    cd bin
    wget https://cnsgenomics.com/software/gcta/bin/gcta_1.93.2beta.zip
    unzip gcta_1.93.2beta.zip && rm gcta_1.93.2beta.zip
    if [ -d gcta ]; then
        mv gcta_1.93.2beta/* ./gcta
    else
        mv gcta_1.93.2beta gcta
    fi
    )
fi

# Install GCTB
# TODO: add links to the binaries, or move binaries to ./bin
if [ ! -f ./bin/gctb/gctb ]; then
    >&2 echo "Downloading GCTB 2.0.3 binaries"
    >&2 echo "Warning: For MPL support, GCTB has to be compiled from source!"
    (
    cd bin
    wget https://cnsgenomics.com/software/gctb/download/gctb_2.03beta_Linux.zip
    unzip gctb_2.03beta_Linux.zip && rm gctb_2.03beta_Linux.zip
    if [ -d gctb ]; then 
        mv gctb_2.03beta_Linux/* gctb
    else
        mv gctb_2.03beta_Linux gctb
    fi
    )
fi

# Download lassosum
if [ ! -d ./bin/lassosum ]; then
    >&2 echo "Downloading Lassosum"
    (
    mkdir bin/lassosum
    cd bin/lassosum
    wget https://github.com/tshmak/lassosum/releases/download/v0.4.5/lassosum_0.4.5.tar.gz
    )
fi

# "Install" GenoPred
if [ ! -d ./workflow/scripts/GenoPred ]; then
    >&2 echo "Downloading GenoPred"
    (
    cd workflow/scripts
    # git clone https://github.com/opain/GenoPred.git && cd GenoPred  && git checkout ${GENOPRED_VERSION}
    git clone git@github.com:intervene-EU-H2020/GenoPred.git && cd GenoPred
    if [ "${GENOPRED_VERSION}" = "latest" ]; then
        git pull
    else
        git checkout ${GENOPRED_VERSION}
    fi
    )
else
    (
    cd ./workflow/scripts/GenoPred
    if [ "${GENOPRED_VERSION}" = "latest" ]; then
         git pull
    else
         git checkout ${GENOPRED_VERSION}
    fi
    )
fi 

# "Install" LDSC
if [ ! -d ./workflow/scripts/ldsc ]; then
    >&2 echo "Downloading LDSC"
    (
    cd workflow/scripts
    git clone https://github.com/bulik/ldsc.git && cd ldsc  && git checkout ${LDSC_VERSION}
    )
else
   (
   cd ./workflow/scripts/ldsc && git checkout ${LDSC_VERSION}
   )
fi 


# "Install" PRSice-2
if [ ! -f ./bin/PRSice_linux ]; then
    >&2 echo "Downloading PRSice-2" 
    (
    cd bin
    wget https://github.com/choishingwan/PRSice/releases/download/${PRSICE_VERSION}/PRSice_linux.zip && unzip PRSice_linux.zip && rm TOY_BASE_GWAS.assoc TOY_TARGET_DATA.bed TOY_TARGET_DATA.bim TOY_TARGET_DATA.fam TOY_TARGET_DATA.pheno
    )
fi
