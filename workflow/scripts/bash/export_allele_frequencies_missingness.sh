#!/bin/bash

shopt -s expand_aliases

# this script is meant to be run in the main analysis folder after harmonization of the target genotypes and ancestry scoring.
# it will create a .tar.gz file that you can then share

if ! command -v plink2 &> /dev/null; then
    echo "could not find plink2 command, trying bin/plink2"
    alias plink2=bin/plink2
fi

set -e

if [[ -z "${BIOBANK}" ]]; then
    echo "Error: BIOBANK undefined, please set the BIOBANK environment variable, e.g.: export BIOBANK='ukbb'"
    exit 1
fi

keepfiles=results/${BIOBANK}/Ancestry_identifier/outlier_detection/*keep

echo "found keep files: ${keepfiles}"

mkdir -p results/${BIOBANK}/afreq

for k in $keepfiles; do

    if [[ $(wc -l <results/ukbb/Ancestry_identifier/outlier_detection/AllAncestry.QC.EUR.keep) -ge 20 ]]; then

        outfile=$(basename ${k})
        outfile=results/${BIOBANK}/afreq/${outfile%%.keep}

        echo "Using output prefix '${outfile}'"

        for chrom in {1..22}; do
            if [ ! -f ${outfile}_chr${chrom}.acount ]; then
                plink2 --bfile custom_input/ukbb/genotypes/chr${chrom} --keep ${k} --missing variant-only --freq counts --out ${outfile}_chr${chrom} --memory 7500
            else
                echo "output file already exists: ${outfile}_chr${chrom}.acount, skipping."
            fi
        done
    else
        echo "fewer than 20 individuals inside {$keepfile}, skipping."
    fi
done

echo "packing..."

tar -zcvf ${BIOBANK}_hapmap3_1kg_allele_frequencies.tar.gz results/${BIOBANK}/afreq

echo "Done. Allele frequency files are packaged in ${BIOBANK}_hapmap3_1kg_allele_frequencies.tar.gz"

