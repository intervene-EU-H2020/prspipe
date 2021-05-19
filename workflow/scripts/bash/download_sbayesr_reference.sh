#!/bin/bash

set -e

# In the future, this script could possibly download reference data for all the different superpopulations
# currently, it only supports the EUR superpopulation provided for GCTB
if [ $# -ne 2 ]; then
    echo "Error: need two input arguments. Superpopulation (e.g. EUR) and output directory."
    exit 1
fi

pop=$1
outdir=$2

if [[ "${pop}" != "EUR" ]]; then
    echo "Error: currently only EUR superpopulation is supported!"
    exit 1
fi

cd ${outdir}
wget -O ukbEURu_hm3_sparse.zip --retry-on-http-error=429 --timeout=5 --waitretry=25 https://zenodo.org/record/3350914/files/ukbEURu_hm3_sparse.zip?download=1


if [ -d ukbEURu_hm3_shrunk_sparse ]; then
    rm -r ukbEURu_hm3_shrunk_sparse
fi

unzip ukbEURu_hm3_sparse.zip && mv ukbEURu_hm3_shrunk_sparse/* ./ && rm -r ukbEURu_hm3_shrunk_sparse/ && rm ukbEURu_hm3_sparse.zip

# Rename matrices so chromosome number is at the end of the file name
# currently these are the same names as provided in GenoPred (step 4.6, https:/opain.github.io/GenoPred/Pipeline_prep.html#46_Prepare_score_and_scale_files_for_polygenic_scoring_using_SBayesR), could be renamed at some point...
for chr in $(seq 1 22);do
mv ukbEURu_hm3_chr${chr}_v3_50k.ldm.sparse.bin ukbEURu_hm3_v3_50k_chr${chr}.ldm.sparse.bin
mv ukbEURu_hm3_chr${chr}_v3_50k.ldm.sparse.info ukbEURu_hm3_v3_50k_chr${chr}.ldm.sparse.info
mv ukbEURu_hm3_chr${chr}_v3_50k_sparse.log ukbEURu_hm3_v3_50k_sparse_chr${chr}.log
done

