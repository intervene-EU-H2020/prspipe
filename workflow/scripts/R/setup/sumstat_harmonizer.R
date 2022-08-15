#!/usr/bin/Rscript

# Checks if the SNP names in the given summary statistics file are in
# RSid format, and does the mapping if they are not in the right format
#
# use --force to force reformatting even if rsids are present (not necessary, but can increase overlap a little)
#
# Author: Sophie Wharrie
# adapted by Remo Monti

suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics file [required]"),
  make_option("--map", action="store", default=NA, type='character',
              help="Path to map of SNP names [required]"),
  make_option("--gz", action="store", default=T, type='logical',
              help="Set to T to gzip summary statistics [optional]"),
  make_option("--output", action="store", default=NA, type='character',
              help="Path for output files [required]"),
  make_option("--force", action="store", type="logical", default=F),
  make_option("--tmpdir", action="store", default=NA, type="character")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

tempdir <- opt$tmpdir
if( is.na(tempdir) ){
    tempdir <- tempdir()
} else {
    if(!dir.exists(tempdir)){
    stop(paste0('specified temporary directory "', tempdir,'" does not exist.'))
    }
}

# open the summary statistics file
if(substr(opt$sumstats,(nchar(opt$sumstats)+1)-3,nchar(opt$sumstats)) == '.gz'){
  GWAS<-data.frame(fread(cmd=paste0('zcat ',opt$sumstats), tmpdir=tempdir))
} else {
  GWAS<-data.frame(fread(opt$sumstats, tmpdir=tempdir))
}

col_names=names(GWAS)

print(paste0("Found ", nrow(GWAS), " SNPs"))

# check the summary statistics file has the expected columns
if(!('SNP' %in% colnames(GWAS))){
  stop("Summary statistics file is missing the SNP column")
}
if(!('POS' %in% colnames(GWAS))){
  stop("Summary statistics file is missing the POS column")
}

if(nrow(GWAS[!(GWAS$SNP %like% "rs"), ])>0){
  print("SNP names not in RSid format... running conversion")
    covert <- T
} else {
  print("SNP names already in RSid format.")
}

if (opt$force){
  print("--force specified, forcing renaming and mapping of SNPs to hg19")
  convert <- T
}


if(convert){
   # open the map file
  if(substr(opt$map,(nchar(opt$map)+1)-3,nchar(opt$map)) == '.gz'){
    map<-data.frame(fread(cmd=paste0('zcat ',opt$map), tmpdir=tempdir))
  } else {
    map<-data.frame(fread(opt$map, tmpdir=tempdir))
  }
  
  # add IUPAC column
  GWAS$IUPAC <- ""
  GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='A']<-'W'
  GWAS$IUPAC[GWAS$A1 == 'C' & GWAS$A2 =='G' | GWAS$A1 == 'G' & GWAS$A2 =='C']<-'S'
  GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='G' | GWAS$A1 == 'G' & GWAS$A2 =='A']<-'R'
  GWAS$IUPAC[GWAS$A1 == 'C' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='C']<-'Y'
  GWAS$IUPAC[GWAS$A1 == 'G' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='G']<-'K'
  GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='C' | GWAS$A1 == 'C' & GWAS$A2 =='A']<-'M'
  GWAS<-GWAS[(GWAS$IUPAC %in% c('R', 'Y', 'K', 'M')),]
  
  # check if there is more overlap for hg38 or hg19
  overlap_hg38 <- nrow(merge(GWAS,map,by.x=c("POS","IUPAC"),by.y=c("pos_hg38","IUPAC_1kg")))
  overlap_hg19 <- nrow(merge(GWAS,map,by.x=c("POS","IUPAC"),by.y=c("pos_hg19","IUPAC_1kg")))

  stopifnot((overlap_hg38 > 0) | (overlap_hg19 > 0))
    
  print(paste0('hg19 overlap: ', overlap_hg19))
  print(paste0('hg38 overlap: ', overlap_hg38))
    
  is_hg38 <- overlap_hg38 > overlap_hg19
    
  if(is_hg38){
    GWAS<-merge(GWAS,map[c("pos_hg38","IUPAC_1kg","pos_hg19","rsid","chr")],by.x=c("POS","IUPAC","CHR"),by.y=c("pos_hg38","IUPAC_1kg","chr"))
    GWAS<-subset(GWAS, select = -c(POS,IUPAC,SNP) )
    names(GWAS)[names(GWAS) == 'pos_hg19'] <- 'POS' # update POS to use hg19
    names(GWAS)[names(GWAS) == 'rsid'] <- 'SNP'
  }
  else{
    GWAS<-merge(GWAS,map[c("IUPAC_1kg","pos_hg19","rsid","chr")],by.x=c("POS","IUPAC","CHR"),by.y=c("pos_hg19","IUPAC_1kg","chr"))
    GWAS<-subset(GWAS, select = -c(IUPAC,SNP) )
    names(GWAS)[names(GWAS) == 'rsid'] <- 'SNP'
  }
  
  GWAS<-GWAS[, col_names]   
}
    
print(paste0("Keeping ", nrow(GWAS), " SNPs"))

# save the modified GWAS file
fwrite(GWAS, opt$output, sep='\t')

if(opt$gz == T){
  # compress the GWAS file
  system(paste0('gzip ', opt$output))
}
