#!/usr/bin/env Rscript

# original author: Oliver Pain
# adapted by: Remo Monti

# Read in environmental variables
source('workflow/scripts/R/source_config.R')


args = commandArgs(trailingOnly=TRUE)
# first argument: superpopulation
# second argument: chunksize
# third argument: chromosome
popul = args[1]
chunksize = as.numeric(args[2])
chrom=as.numeric(args[3])

options(scipen=999)

# Create list of SNPs that have MAF above 0.0 in the EUR 1KG sample
for(i in chrom){
  frq<-read.table(paste0(Geno_1KG_dir, '/freq_files/',popul,'/1KGPhase3.w_hm3.',popul,'.chr',i,'.frq'), header=T)
  frq<-frq[frq$MAF > 0.0,]
  write.table(frq$SNP,paste0('resources/LD_matrix/sbayesr/1000G/fromscratch/',popul,'/chr',i,'/SNP.txt'), col.names=F, row.names=F, quote=F)
  nsnp <- nrow(frq)
  if (nsnp < chunksize){
    nsnp_chunk <- 1
  } else {
    nsnp_chunk<-ceiling(nsnp/chunksize)
  }
  for(j in 1:nsnp_chunk){
    start<-(chunksize*(j-1))+1
    end<-chunksize*j
    # print(start)
    # print(end)
    outfile <- paste0('resources/LD_matrix/sbayesr/1000G/fromscratch/',popul,'/chr',i,'/SNP_snp',start,'-',end,'.txt')
    sink(outfile)
    cat(start,' ',end,'\n')
    sink()
  }
}