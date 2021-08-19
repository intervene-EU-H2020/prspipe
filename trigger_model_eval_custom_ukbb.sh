#!/bin/bash


if [ $# -gt 2 ]; then
    echo "Error: script only takes two arguments : Phenotype [-n|--dryrun]"
    exit 1
fi


DRYRUN=""
if [ $# -gt 1 ]; then
    if [ "$2" != "-n" ]; then
        if [ "$2" != "--dryrun" ]; then
            echo "Error: second argument can only be [-n|--dryrun]"
            echo "found: ${2}"
            exit 1
        else
            DRYRUN='-n'
        fi
    else
        DRYRUN='-n'
    fi
fi

GRP="$(ls -1 custom_input/predictor_groups/*.predictor_groups | cut -d/ -f3 | sed 's/\.predictor_groups//g')"

OUTFILES="$(echo ${GRP} | tr ' ' '\n' | awk -v pheno=$1 '{print "./results/UKBB/evaluation/" $1 "/" pheno "/best_models.tsv" }' | tr '\n' ' ')"

bash run.sh ${DRYRUN} ${OUTFILES} 