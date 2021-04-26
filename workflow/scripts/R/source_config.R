#!/usr/bin/env Rscript

source_config <- function(){
    config <- snakemake@config    
    for (i in seq_along(config)){
        assign(names(config)[i],config[[i]],envir = .GlobalEnv)
    }
}

source_config()

print(plink1_9)


