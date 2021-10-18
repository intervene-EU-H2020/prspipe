#!/usr/bin/env Rscript

library(data.table)

args = commandArgs(trailingOnly=TRUE)

stopifnot(length(args) == 6)

# hg18 bim
hg18_bim = args[1]
# output of liftOver
hg19_lo = args[2]
hg38_lo = args[3]
# the chromosome
chrom = args[4]
# the dbSNP version
dbsnp_version = args[5]
# the output file
outfile = args[6]

# load dbSNP from the custom location
.libPaths(c(.libPaths(), "R/local/"))
library('GenomicRanges')
library(paste0("SNPlocs.Hsapiens.dbSNP",dbsnp_version,".GRCh38"),character.only = TRUE)

matchtable <- fread(hg18_bim, sep = '\t', col.names = c('chr_hg18','rsid','dst','pos_hg18','a1','a2'))

# merge with hg19
matchtable <- merge(matchtable, fread(hg19_lo, sep=' '), by = c('rsid', 'a1', 'a2'))
setnames(matchtable,'chr','chr_hg19')

# merge with hg38
matchtable <- merge(matchtable, fread(hg38_lo, sep=' '), by = c('rsid', 'a1', 'a2'), suffixes = c('','_hg38'))
setnames(matchtable, 'chr','chr_hg38')

# discard ambiguous chromosome SNPs
matchtable <- matchtable[chr_hg19 == chr_hg38]
matchtable <- matchtable[chr_hg18 == chr_hg38]

# sort
matchtable <- matchtable[order(chr_hg19, pos_hg19, a1, a2)]

# remove redundant columns
matchtable[,chr_hg19:=NULL]
matchtable[,chr_hg38:=NULL]
setnames(matchtable, 'chr_hg18','chr')

# filter by chromosome
matchtable <- matchtable[chr == chrom]

# reorder columns
setcolorder(matchtable, c('chr', 'rsid', 'pos_hg18', 'pos_hg19', 'pos_hg38', 'a1', 'a2'))

# add IUPAC codes
matchtable[,IUPAC:='']
matchtable[ (a1 == 'A' & a2 =='T') | (a1 == 'T' & a2 =='A'), IUPAC:='W']
matchtable[ (a1 == 'C' & a2 =='G') | (a1 == 'G' & a2 =='C'), IUPAC:='S']
matchtable[ (a1 == 'A' & a2 =='G') | (a1 == 'G' & a2 =='A'), IUPAC:='R']
matchtable[ (a1 == 'C' & a2 =='T') | (a1 == 'T' & a2 =='C'), IUPAC:='Y']
matchtable[ (a1 == 'G' & a2 =='T') | (a1 == 'T' & a2 =='G'), IUPAC:='K']
matchtable[ (a1 == 'A' & a2 =='C') | (a1 == 'C' & a2 =='A'), IUPAC:='M']

options(scipen = 999)

snp_pos <- GRanges(seqnames=matchtable$chr, ranges=IRanges(start=matchtable$pos_hg38, width=1), rsid=matchtable$rsid)

cat(paste0('Overlapping ',length(snp_pos),' SNP positions with dbSNP annotation ',dbsnp_version,'.\n'))

# get the overlaps by position with dbSNP
snp_overlaps <- snpsByOverlaps(get(paste0('SNPlocs.Hsapiens.dbSNP',dbsnp_version,'.GRCh38')), ranges = snp_pos)
cat(paste0('Found ', length(snp_overlaps), ' overlapping SNPs.\n'))
    
snp_overlaps<-data.table(as.data.frame(snp_overlaps))
setnames(snp_overlaps, old=c('seqnames','pos','strand','RefSNP_id','alleles_as_ambig'), new=c('chr','pos_hg38','strand','rsid','IUPAC'))
snp_overlaps[,strand:=NULL]
matchtable$chr <- as.character(matchtable$chr)
snp_overlaps$chr <- as.character(snp_overlaps$chr)

matchtable <- merge(matchtable, snp_overlaps, by = c('chr','pos_hg38','IUPAC'), suffixes = c('_hm3','_dbSNP'), all.x = TRUE)
stopifnot(nrow(matchtable) > 0) # the intersection is empty. Something is wrong...

# get rid of ambiguous SNPs, and replace IDs with latest versions...
matchtable[,keep:=FALSE]
matchtable[,rsid_mrg:='.']

# duplicated, but at least one ID matches with dbSNP -> keep the one with matching ID
matchtable[(duplicated(matchtable[,list(chr, pos_hg38, IUPAC)], fromLast = T) | duplicated(matchtable[,list(chr, pos_hg38, IUPAC)], fromLast = F)) & (rsid_hm3 == rsid_dbSNP), rsid_mrg:=rsid_dbSNP]
matchtable[(duplicated(matchtable[,list(chr, pos_hg38, IUPAC)], fromLast = T) | duplicated(matchtable[,list(chr, pos_hg38, IUPAC)], fromLast = F)) & (rsid_hm3 == rsid_dbSNP), keep:=TRUE]

# not duplicated -> assign the new ID
matchtable[!(duplicated(matchtable[,list(chr, pos_hg38, IUPAC)], fromLast = T) | duplicated(matchtable[,list(chr, pos_hg38, IUPAC)], fromLast = F)) & (!is.na(rsid_dbSNP)), rsid_mrg:=rsid_dbSNP]
matchtable[!(duplicated(matchtable[,list(chr, pos_hg38, IUPAC)], fromLast = T) | duplicated(matchtable[,list(chr, pos_hg38, IUPAC)], fromLast = F)) & (!is.na(rsid_dbSNP)), keep:=TRUE]

# ID not listed in dbSNP (see restriction of variants in the R-package...) -> assign ID given in HapMap3 (might be old...)
matchtable[is.na(rsid_dbSNP),rsid_mrg:=rsid_hm3]
matchtable[is.na(rsid_dbSNP),keep:=TRUE]

# remove rows not marked with "keep", and drop unnecessary columns
cat(paste0('Removing ',sum(!matchtable$keep),' matches because of conflicts.'))
matchtable <- matchtable[(keep)]
matchtable[,keep:=NULL]
matchtable[,dst:=NULL]
                               
# export
fwrite(matchtable, outfile, sep='\t', scipen=50, na='.', quote=FALSE)