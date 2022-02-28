#!/usr/bin/env Rscript

# Use R to identify a list of SNPs that matches with the hapmap3 snplist

# script derived from the GenoPred pipeline
# https://opain.github.io/GenoPred/
# original author: Oliver Pain
# adapted by Remo Monti

# Read in environmental variables
source('workflow/scripts/R/source_config.R')
Geno_1KG_dir <- 'resources/1kg'

library(data.table)
hapmap3_snps<-fread(paste0(HapMap3_snplist_dir,'/w_hm3.snplist'))

# Generate IUPAC codes
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'A' & hapmap3_snps$A2 =='T' | hapmap3_snps$A1 == 'T' & hapmap3_snps$A2 =='A']<-'W'
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'C' & hapmap3_snps$A2 =='G' | hapmap3_snps$A1 == 'G' & hapmap3_snps$A2 =='C']<-'S'
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'A' & hapmap3_snps$A2 =='G' | hapmap3_snps$A1 == 'G' & hapmap3_snps$A2 =='A']<-'R'
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'C' & hapmap3_snps$A2 =='T' | hapmap3_snps$A1 == 'T' & hapmap3_snps$A2 =='C']<-'Y'
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'G' & hapmap3_snps$A2 =='T' | hapmap3_snps$A1 == 'T' & hapmap3_snps$A2 =='G']<-'K'
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'A' & hapmap3_snps$A2 =='C' | hapmap3_snps$A1 == 'C' & hapmap3_snps$A2 =='A']<-'M'

hapmap3_snps$SNP_IUPAC<-paste(hapmap3_snps$SNP,hapmap3_snps$IUPAC,sep=':')

# For each chr
for(chr in 1:22){
    bim<-fread(paste0(Geno_1KG_dir,'/1KGPhase3.chr',chr,'.bim'))

    bim$IUPAC[bim$V5 == 'A' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='A']<-'W'
    bim$IUPAC[bim$V5 == 'C' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='C']<-'S'
    bim$IUPAC[bim$V5 == 'A' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='A']<-'R'
    bim$IUPAC[bim$V5 == 'C' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='C']<-'Y'
    bim$IUPAC[bim$V5 == 'G' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='G']<-'K'
    bim$IUPAC[bim$V5 == 'A' & bim$V6 =='C' | bim$V5 == 'C' & bim$V6 =='A']<-'M'

    bim$SNP_IUPAC<-paste(bim$V2,bim$IUPAC,sep=':')
    
    bim_hapmap3_snps<-merge(hapmap3_snps,bim,by='SNP_IUPAC')
    sum(duplicated(bim_hapmap3_snps$V2))
    
    bim$V2[!(bim$SNP_IUPAC %in% bim_hapmap3_snps$SNP_IUPAC)] <- paste0(bim$V2[!(bim$SNP_IUPAC %in% bim_hapmap3_snps$SNP_IUPAC)],'_excl')
    bim[!(bim$SNP_IUPAC %in% bim_hapmap3_snps$SNP_IUPAC),]
    
    fwrite(bim[,1:6], paste0(Geno_1KG_dir,'/1KGPhase3.chr',chr,'.bim'), sep='\t', col.names=F)
    fwrite(as.list(bim_hapmap3_snps$V2), paste0(Geno_1KG_dir,'/1KGPhase3.chr',chr,'.extract'), sep='\n', col.names=F)   
}

warnings()
