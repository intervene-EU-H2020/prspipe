#!/usr/bin/env Rscript

# We have already calculated ancestry PCs for UK biobank and identified European ancestry based on 1KG Phase participants.
# Now create a random subset of 10K European individuals

# script derived from the GenoPred pipeline
# https://opain.github.io/GenoPred/
# original author: Oliver Pain
# adapted by Remo Monti

# Note: Unlike in the original GenoPred script we don't care if the MAF/LD is calculated on the target individuals in the UKBB

# Read in environmental variables
source('workflow/scripts/R/source_config.R')

library(data.table)

iid_withdraw <- as.character(read.table(UKBB_exclude, header=FALSE)[,'V1'])

# Read in list of EUR UKBB individuals
EUR_UKBB<-fread(paste0(UKBB_output, '/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.EUR.keep'))
names(EUR_UKBB)<-c('FID','IID')

EUR_UKBB[,FID:=as.character(FID)]
EUR_UKBB[,IID:=as.character(IID)]

EUR_UKBB<-EUR_UKBB[!(IID %in% iid_withdraw)]

# Extract a random 10K individuals
set.seed(1)
EUR_UKBB_10K<-EUR_UKBB[sample(dim(EUR_UKBB)[1],10000),]

write.table(EUR_UKBB_10K, paste0(UKBB_output, '/UKBB_ref/keep_files/UKBB_noPheno_EUR_10K.keep'), col.names=F, row.names=F, quote=F)

