#!/bin/bash

# script that downloads the PRS used for the methods comparison (May 2022)

echo "downloading data for methods comparison, package 1 (may 2022)"


if ! command -v wget &> /dev/null
then
    echo "could not find wget, trying to use curl"

    download () {
        curl -L -o $1 $2
    }
    

else
    echo "using wget"

    download () {
        wget -O $1 $2
    }
    
fi


echo "downloading pre-computed PRS"

if [ ! -f prs.tar.gz ]; then
    download prs.tar.gz 'https://figshare.com/ndownloader/files/35016199?private_link=7f738b939e1cba580708'
else
   echo "skipped: prs.tar.gz already exists."
fi


if [ $? -ne 0 ]; then
    echo "Error while trying to download prs!"
    rm -f prs.tar.gz
    exit 1
fi


echo "downloading liftOver"

mkdir -p liftover

(
cd liftover
if [ ! -f liftOver ]; then
    download liftOver https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && \
    chmod u+x liftOver && \
    download hg18ToHg19.over.chain.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz && \
    gunzip hg18ToHg19.over.chain.gz && \
    download hg18ToHg38.over.chain.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz && \
    gunzip hg18ToHg38.over.chain.gz && \
    download hg19ToHg38.over.chain.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz && \
    gunzip hg19ToHg38.over.chain.gz 
else
    echo "skipped: liftover/liftOver already exists."
fi
)

if [ $? -ne 0 ]; then
    echo "Warning: Errors encountered when trying to download liftOver"
fi

echo "downloading 1000 genomes metadata"

mkdir -p resources/1kg

(
if [ ! -f resources/1kg/integrated_call_samples_v3.20130502.ALL.panel_small ]; then
    cd resources/1kg && \
    download integrated_call_samples_v3.20130502.ALL.panel ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel && \
    cut -f 1-3 integrated_call_samples_v3.20130502.ALL.panel > integrated_call_samples_v3.20130502.ALL.panel_small
else
    echo "skipped: resources/1kg/integrated_call_samples_v3.20130502.ALL.panel_small already exists."
fi
)

if [ $? -ne 0 ]; then
    echo "Error when trying to download 1000 genomes metadata!"
    rm -f resources/1kg/integrated_call_samples_v3.20130502.ALL.panel_small
    exit 1
fi

echo "downloading pre-processed 1000 genomes reference"

if [ ! -f 1KGPhase3.w_hm3.chr.tar.gz ]; then
    download 1KGPhase3.w_hm3.chr.tar.gz 'https://figshare.com/ndownloader/files/37059184?private_link=cea777bb772ed1dc4ca8'
else
    echo "skipped: 1KGPhase3.w_hm3.chr.tar.gz already exists."
fi

if [ $? -ne 0 ]; then
    echo "Error when downloading pre-processed 1000 genomes reference!"
    rm 1KGPhase3.w_hm3.chr.tar.gz
fi

echo "unpacking pre-processed 1000 genomes reference"

tar -xzvf 1KGPhase3.w_hm3.chr.tar.gz && rm 1KGPhase3.w_hm3.chr.tar.gz && \
find resources/1kg/ -name '1KGPhase3.w_hm3.chr*' -type f -exec touch {} + 

if [ $? -ne 0 ]; then
    echo "Error unpacking pre-processed 1000 genomes reference!"
    exit 1
fi

