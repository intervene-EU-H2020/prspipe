#!/bin/bash

# Date:         24.03.2021
# Author:       Remo Monti & Sophie Wharrie
# Description:  Download binaries, move them to bin, clone GenoPred repo

GENOPRED_VERION="4a0eaec89acb336c532a59dd005e914704184d32"
LDSC_VERSION="aa33296abac9569a6422ee6ba7eb4b902422cc74"

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
if [ ! -f ./bin/lassosum ]; then
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
    git clone https://github.com/opain/GenoPred.git && cd GenoPred  && git checkout ${GENOPRED_VERSION}
    )
else
   cd ./workflow/scripts/GenoPred && git checkout ${GENOPRED_VERSION}
fi 

# "Install" LDSC
if [ ! -d ./workflow/scripts/ldsc ]; then
    >&2 echo "Downloading LDSC"
    (
    cd workflow/scripts
    git clone https://github.com/bulik/ldsc.git && cd ldsc  && git checkout ${LDSC_VERSION}
    )
else
   cd ./workflow/scripts/ldsc && git checkout ${LDSC_VERSION}
fi 



