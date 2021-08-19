#!/bin/bash

# 1st argument: comma-separated list of study ids
# 2nd argument: comma-separated list if prs methods (options: ldpred, ldpred2, lassosum, prscs, sblup, sbayesr, dbslmm, pt_clump)
# 3rd argument: name of predictor group

# writes a predictor_groups-file to custom_input/predictor_groups/ which can be used to run Model_Builder_V2.R

if [ $# -lt 3 ]; then
    echo "usage: write_predictor_groups_ukbb.sh studies methods groupname [pheno] [-y] "
    exit 1
fi 

if [ $# -eq 5 ]; then
    if [[ "$5" != '-y' ]]; then
        echo "fifth argument can only be '-y' ( yes -> run snakemake)"
        exit 1
    fi
fi


STUDIES=$1
METHODS=$2
NAME=$3

outfile="custom_input/predictor_groups/${NAME}.predictor_groups"

if [ -f $outfile ]; then
    echo "Error: output file ${outfile} already exists!"
    exit 1
else
    echo "writing predictors to ${outfile}."
fi

METHODS="$(echo $METHODS | tr ',' '|' )"

echo 'predictors group' >> ${outfile}

for study in $(echo $STUDIES | tr ',' '\n'); do
    egrep $METHODS ./results/UKBB/PRS_for_comparison/evaluation/${study}/Association_withPRS/UKBB.w_hm3.${study}.AllMethodComp.predictor_groups | awk -v s=${study} '{print $1, $2"."s}' >> ${outfile}
done

if [ $# -ge 4 ]; then
    PHENO=$4
    if [[ $5 == '-y' ]]; then
        ./run.sh -n "./results/UKBB/evaluation/${NAME}/${PHENO}/AllMethodComp.ok"
    else
        echo "./results/UKBB/evaluation/${NAME}/${PHENO}/AllMethodComp.ok"
    fi
fi

