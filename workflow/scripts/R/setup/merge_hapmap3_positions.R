#!/usr/bin/env Rscript

library(data.table)

args = commandArgs(trailingOnly=TRUE)

stopifnot(length(args) == 4)

# hg18 bim
hg18_bim = args[1]
# output of liftOver
hg19_lo = args[2]
hg38_lo = args[3]
# the output file
outfile = args[4]

matchtable <- fread(hg18_bim, sep = '\t', col.names = c('chr_hg18','rsid','dist','pos_hg18','a1','a2'))

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

# export
fwrite(matchtable, outfile, sep='\t', scipen=50)