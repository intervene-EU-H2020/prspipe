#!/usr/bin/env Rscript

# We need to divide the reference genetic data into LD blocks. Previously defined LD blocks have been downloaded from here: https://bitbucket.org/nygcresearch/ldetect-data/src/master/

# First, we need to create a SNPlist for each LD block.

# script derived from the GenoPred pipeline
# https://opain.github.io/GenoPred/
# original author: Oliver Pain
# adapted by Remo Monti

# Read in environmental variables
source('workflow/scripts/R/source_config.R')

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

# 1000 genomes superpopulation
popul <- args[1]

# Read in bim file for reference
bim_all<-NULL
for(i in 1:22){
    bim<-fread(paste0(Geno_1KG_dir,'/1KGPhase3.w_hm3.chr',i,'.bim'))
    bim_all<-rbind(bim_all, bim)
}

# Read in the freq files to create the snpinfo file
frq_all<-NULL
for(i in 1:22){
    frq<-fread(paste0(Geno_1KG_dir,'/freq_files/',popul,'/1KGPhase3.w_hm3.',popul,'.chr',i,'.frq'))
    frq_all<-rbind(frq_all, frq)
}

# Create snpinfo file
frq_bim<-merge(bim_all, frq_all[,c('SNP','MAF','NCHROBS')], by.x='V2', by.y='SNP')
names(frq_bim)<-c('SNP','CHR','POS','BP','A1','A2','MAF','NCHROBS')
frq_bim<-frq_bim[,c('CHR','SNP','BP','A1','A2','MAF','NCHROBS')]
frq_bim<-frq_bim[order(frq_bim$CHR, frq_bim$BP),]

# Remove SNPs with a MAF < 0.01
frq_bim<-frq_bim[frq_bim$MAF > 0.01,]
# Remove SNPs with a missingness > 0.01
frq_bim<-frq_bim[frq_bim$NCHROBS > (max(frq_bim$NCHROBS)*0.99),]
frq_bim$NCHROBS<-NULL

# Extract snplist for each range in LD Block bed files
bed<-fread('resources/ldetect-data/',popul,'/fourier_ls-all.bed')
bed$chr<-as.numeric(gsub('chr','',bed$chr))

blk_chr<-bed$chr
write.table(blk_chr,paste0('/resources/LD_matrix/prscs/1000G/fromscratch/',popul,'/blk_chr'), col.names=F, row.names=F, quote=F)

blk_size<-NULL
snpinfo_1kg_hm3<-NULL
for(i in 1:dim(bed)[1]){
  frq_bim_i<-frq_bim[frq_bim$CHR == bed$chr[i] & frq_bim$BP > bed$start[i] & frq_bim$BP < bed$stop[i],]
  blk_size<-c(blk_size,dim(frq_bim_i)[1])
  snpinfo_1kg_hm3<-rbind(snpinfo_1kg_hm3,frq_bim_i)
  write.table(frq_bim_i$SNP,paste0('/resources/LD_matrix/prscs/1000G/fromscratch/',popul,'/LD_Blocks/Block_',i,'.snplist'), col.names=F, row.names=F, quote=F)
}

write.table(blk_size,paste0('/resources/LD_matrix/prscs/1000G/fromscratch/',popul,'/LD_Blocks/blk_size'), col.names=F, row.names=F, quote=F)

write.table(snpinfo_1kg_hm3, paste0('/resources/LD_matrix/prscs/1000G/fromscratch/',popul,'/LD_blocks/snpinfo_1kg_hm3'), col.names=T, row.names=F, quote=F)
