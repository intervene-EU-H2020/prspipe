#!/usr/bin/env Rscript

library(bigsnpr)
library(bigreadr)

# GenoPred 4.8.1 Create LD reference for LDPred2
# Here we create reference files for polygenic scores calculated by LDPred2, a method for performing bayesian shrinkage analysis with summary data and an LD-reference
# original author: Oliver Pain
# adapted by: Remo Monti

# read configuration
source('workflow/scripts/R/source_config.R')

args = commandArgs(trailingOnly=TRUE)

# first argument: population
popul = args[1]
# second argument: number of cores
if (length(args)>1){
  NCORES = as.integer(args[2])
  } else {
  NCORES = 6
}


# Subset Europeans from reference
system(paste0(plink1_9,' --bfile ',Geno_1KG_dir,'/1KGPhase3.w_hm3.GW --keep ',Geno_1KG_dir,'/keep_files/',popul,'_samples.keep --make-bed --out resources/LD_matrix/ldpred2/1000G/fromscratch/',popul,'/1KG.',popul,'.tmp.GW'))

# Read in reference data
snp_readBed(paste0('resources/LD_matrix/ldpred2/1000G/fromscratch/',popul,'/1KG.',popul,'.tmp.GW.bed'))

# Attach the ref object in R session
ref <- snp_attach(paste0('resources/LD_matrix/ldpred2/1000G/fromscratch/',popul,'/1KG.',popul,'.tmp.GW.rds'))
G <- ref$genotypes
bigassertr::assert_dir(paste0('resources/LD_matrix/ldpred2/1000G/fromscratch/',popul))

#### Impute missing values (bigsnpr can't handle missing data in most functions)
G_imp<-snp_fastImputeSimple(G, method = "mean2", ncores = NCORES)
G<-G_imp

# Save imputed reference
ref$genotypes<-G
saveRDS(ref, paste0(Geno_1KG_dir,'/1KGPhase3.w_hm3.GW.rds'))

#### Compute LD matrices ####
CHR <- ref$map$chr
POS <- ref$map$physical.pos
POS2 <- snp_asGeneticPos(CHR, POS, dir =paste0(genetic_map,'/CEU'), ncores = NCORES)

# Compute LD
for(chr in 1:22){
  print(chr)
  ind.chr <- which(CHR == chr)
  corr <- snp_cor(G, ind.col = ind.chr, infos.pos = POS2[ind.chr], size = 3 / 1000, ncores = NCORES)
  saveRDS(corr, file = paste0('resources/LD_matrix/ldpred2/1000G/fromscratch/',popul,'/LD_chr', chr, ".rds"), version = 2)
}

# Compute LD scores
ref$map$ld <- do.call('c', lapply(1:22, function(chr) {
  cat(chr, ".. ", sep = "")
  corr_chr <- readRDS(paste0('resources/LD_matrix/ldpred2/1000G/fromscratch/',popul,'/LD_chr', chr, ".rds"))
  Matrix::colSums(corr_chr^2)
}))

saveRDS(ref$map, paste0('resources/LD_matrix/ldpred2/1000G/fromscratch/',popul,'/map.rds'), version = 2)

# Save reference SD of genotypes
sd <- runonce::save_run(
  sqrt(big_colstats(G, ncores = NCORES)$var),
  file = paste0('resources/LD_matrix/ldpred2/1000G/fromscratch/',popul,'/sd.rds')
)