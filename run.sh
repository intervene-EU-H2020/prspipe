#!/bin/bash

snakemake --cores 1 --snakefile workflow/Snakefile "${@}"
