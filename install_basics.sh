#!/bin/bash

# Description:  Download binaries, move them to bin, and clone GenoPred repo

GENOPRED_VERION="latest"
LDSC_VERSION="aa33296abac9569a6422ee6ba7eb4b902422cc74"

if [ ! -f README.md ]; then
    echo "Error: Are you in the project base directory? Abort."
    exit 1
fi

if [ ! -d ./bin ]; then
    mkdir bin
fi


# Download Plink 2.0 
if [ ! -f ./bin/plink2 ]; then
    >&2 echo "Downloading Plink 2.0 binaries"
    (
    cd bin
    wget -O plink2.zip 'https://figshare.com/ndownloader/files/33895304?private_link=3641ef6df51eddbeea60'
    unzip plink2.zip && rm plink2.zip
    avx_support="$(grep avx /proc/cpuinfo | wc -l)"
    # if this fails, check https://www.cog-genomics.org/plink/2.0/ for the latest download links
    if [ $avx_support -gt 0 ]; then
        echo 'found AVX support, will use the plink 2.0 AVX version'
    else
        echo 'no AVX support found, will use the plink 2.0 64-bit version'
        echo 'you can still switch to the avx version later...'
        mv plink2 plink2_avx && mv plink2_x86_64 plink2
    fi
    )
fi

# Install Plink 1.9
if [ ! -f ./bin/plink ]; then
    >&2 echo "Downloading Plink 1.9 binaries"
    (
    cd bin
    # if this fails, check https://www.cog-genomics.org/plink/1.9/ for the latest download links
    wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip
    unzip plink_linux_x86_64_20201019.zip && rm plink_linux_x86_64_20201019.zip
    )
fi


# Install QCTOOL 2
if [ ! -f ./bin/qctool ]; then
    >&2 echo "Downloading QCTOOL V2" 
    (
    # TODO: configure version (?)
    cd bin
    wget https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.6-Ubuntu16.04-x86_64.tgz
    tar -xvf qctool_v2.0.6-Ubuntu16.04-x86_64.tgz && rm qctool_v2.0.6-Ubuntu16.04-x86_64.tgz
    mv qctool_v2.0.6-Ubuntu16.04-x86_64/qctool ./ && rm -r qctool_v2.0.6-Ubuntu16.04-x86_64
    )
fi

# "Install" LDSC
if [ ! -d ./workflow/scripts/ldsc/ ]; then
    >&2 echo "Downloading LDSC"
    (
    cd workflow/scripts
    git clone https://github.com/bulik/ldsc.git && cd ldsc && git checkout ${LDSC_VERSION}
    )
else
   (
   cd ./workflow/scripts/ldsc && git checkout ${LDSC_VERSION}
   )
fi 

# clone GenoPred
if [ ! -d ./workflow/scripts/GenoPred ]; then
    >&2 echo "Downloading GenoPred"
    (
    cd workflow/scripts
    git clone git@github.com:intervene-EU-H2020/GenoPred.git || git clone https://github.com/intervene-EU-H2020/GenoPred.git
    cd GenoPred
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
