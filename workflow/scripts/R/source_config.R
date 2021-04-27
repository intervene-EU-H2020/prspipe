#!/usr/bin/env Rscript

source_config <- function(){
    # load variables defined in config.yaml to global namespace
    if (exists('snakemake')){
        config <- snakemake@config
    } else {
        require(yaml)
        config <- yaml.load_file('config/config.yaml')
    }
    for (i in seq_along(config)){
        assign(names(config)[i],config[[i]],envir = .GlobalEnv)
    }
}

source_config()
