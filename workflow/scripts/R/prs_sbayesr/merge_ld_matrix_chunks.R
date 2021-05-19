#!/usr/bin/env Rscript

# Read in environmental variables
# original author: Oliver Pain
# adapted by: Remo Monti

source('workflow/scripts/R/source_config.R')

args = commandArgs(trailingOnly=TRUE)
popul = args[1] # population
args = args[-1]
chrom = args[1] # chromosome
files = args[-1] # remaining arguments are the .bin files to merge

options(scipen=999)


files<-gsub('.bin', '', files)
file_num<-gsub('.*.snp', '', files)
file_num<-as.numeric(gsub('-.*', '', file_num))
# Sort files in order of genomic location otherwise SBayeR does not converge!
files<-files[order(file_num)]

write.table(files,paste0('resources/LD_matrix/sbayesr/1000G/fromscratch/',popul,'/shrunk_ld_chr',chrom,'merge_list'), col.names=F, row.names=F, quote=F)
if (length(files) == 1){
    # prevents a crash
    system(paste0(gctb,' --ldm ',files[1],' --make-shrunk-ldm --out resources/LD_matrix/sbayesr/1000G/fromscratch/',popul,'/1KGPhase3.w_hm3.',popul,'.chr',chrom))
} else { 
    system(paste0(gctb,' --mldm resources/LD_matrix/sbayesr/1000G/fromscratch/',popul,'/shrunk_ld_chr',chrom,'merge_list --make-shrunk-ldm --out resources/LD_matrix/sbayesr/1000G/fromscratch/',popul,'/1KGPhase3.w_hm3.',popul,'.chr',chrom))
}
