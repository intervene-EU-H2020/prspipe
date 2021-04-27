#!/bin/bash

snakemake --cores 4 --snakefile workflow/Snakefile "${@}"
