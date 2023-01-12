#!/usr/bin/env Rscript

# rule: score_to_target_geno_overlap_create_rsid_rds (workflow/rules/bigsnpr_allele_freq_ancestry.smk)

library(data.table)

score_paths <- fread(commandArgs(trailingOnly=TRUE), header=FALSE)$V1

get_rsid <- function(sspath){
    fread(sspath, select = 'SNP')$SNP
}

var_match <- lapply(score_paths, get_rsid)

names(var_match) <- score_paths

saveRDS(var_match, file = 'temp/prs_rsid.rds')

