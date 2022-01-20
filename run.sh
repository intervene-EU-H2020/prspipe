#!/bin/bash

if [ -z ${SNAKEMAKE_CORES+x} ]; then
    SNAKEMAKE_CORES=1
    # TODO: handle basic configuration like this with profiles https://snakemake.readthedocs.io/en/stable/executing/cli.html?highlight=default%20profile#profiles
fi

snakemake --use-conda --cores "${SNAKEMAKE_CORES}" --snakefile workflow/Snakefile "${@}"

