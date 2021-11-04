#!/usr/bin/env Rscript

suppressMessages(library(data.table))

suppressMessages(library("optparse"))
option_list = list(
    make_option("--mapping", action="store", default=NA, type='character', help="Path to mapping file"),
    make_option("--lifted", action="store", default=NA, type="character", help="Path to lifted file"),
    make_option("--chr", action="store", default=NA, type="character", help="chromosome"),
    make_option("--dbsnp_version_hg37", action="store", default="144", type="character"),
    make_option("--dbsnp_version_hg38", action="store", default="151", type="character"),
    make_option("--out", action="store", default=NA, type='character', help='Output file.')
)

opt = parse_args(OptionParser(option_list=option_list))
# load dbSNP from the custom location
.libPaths(c(.libPaths(), "R/local/"))
suppressMessages(library('GenomicRanges'))
suppressMessages(library(paste0("SNPlocs.Hsapiens.dbSNP",opt$dbsnp_version_hg38,".GRCh38"),character.only = TRUE))


# colnames: chr, pos_hg19, rsid_dbSNP, a1, a2, IUPAC_1kg, IUPAC_dbSNP
matchtable <- fread(opt$mapping, sep = '\t', header=T)
setnames(matchtable,c('chr','rsid_dbSNP','IUPAC_dbSNP'),c('chr_hg19','rsid_dbsnp',paste0('IUPAC_dbSNP', opt$dbsnp_version_hg37)))

# merge with hg38
matchtable <- merge(matchtable, fread(opt$lifted, sep='\t'), by = c('rsid_1kg','rsid_dbsnp', 'a1', 'a2'))
setnames(matchtable,'chr','chr_hg38')

# discard ambiguous chromosome SNPs
print(sum(matchtable$chr_hg19 != matchtable$chr_hg38))
matchtable <- matchtable[chr_hg19 == chr_hg38]

# sort
matchtable <- matchtable[order(chr_hg19, pos_hg19, a1, a2)]

# remove redundant columns
matchtable[,chr_hg38:=NULL]
setnames(matchtable, 'chr_hg19','chr')

# filter by chromosome
matchtable <- matchtable[chr == opt$chr]

# reorder columns
setcolorder(matchtable, c('chr', 'rsid_1kg', 'rsid_dbsnp', 'pos_hg19', 'pos_hg38', 'a1', 'a2'))

options(scipen = 999)

snp_pos <- GRanges(seqnames=matchtable$chr, ranges=IRanges(start=matchtable$pos_hg38, width=1), rsid=matchtable$rsid)

cat(paste0('Overlapping ',length(snp_pos),' SNP positions with dbSNP annotation ',opt$dbsnp_version_hg38,'.\n'))

# get the overlaps by position with dbSNP
snp_overlaps <- snpsByOverlaps(get(paste0('SNPlocs.Hsapiens.dbSNP',opt$dbsnp_version_hg38,'.GRCh38')), ranges = snp_pos)
cat(paste0('Found ', length(snp_overlaps), ' overlapping SNPs.\n'))
    
snp_overlaps<-data.table(as.data.frame(snp_overlaps))
setnames(snp_overlaps, old=c('seqnames','pos','strand','RefSNP_id','alleles_as_ambig'), new=c('chr','pos_hg38','strand','rsid',paste0('IUPAC_dbSNP', opt$dbsnp_version_hg38)))
snp_overlaps[,strand:=NULL]
matchtable$chr <- as.character(matchtable$chr)
snp_overlaps$chr <- as.character(snp_overlaps$chr)

setnames(matchtable, 'rsid_dbsnp', 'rsid')
matchtable <- merge(matchtable, snp_overlaps, by = c('chr','pos_hg38'), suffixes =c('_old','_new'), all.x = TRUE)

stopifnot(nrow(matchtable) > 0) # if the intersection is empty. Something is wrong...

# get rid of ambiguous SNPs, and replace IDs with latest versions...
matchtable[,keep:=FALSE]
matchtable[,rsid_mrg:='.']

# match only position
# duplicated, but at least one ID matches with dbSNP -> keep the one with matching ID
matchtable[(duplicated(matchtable[,list(chr, pos_hg38)], fromLast = T) | duplicated(matchtable[,list(chr, pos_hg38)], fromLast = F)) & (rsid_old == rsid_new), rsid_mrg:=rsid_new]
matchtable[(duplicated(matchtable[,list(chr, pos_hg38)], fromLast = T) | duplicated(matchtable[,list(chr, pos_hg38)], fromLast = F)) & (rsid_old == rsid_new), keep:=TRUE]

# not duplicated -> assign the new ID
matchtable[!(duplicated(matchtable[,list(chr, pos_hg38)], fromLast = T) | duplicated(matchtable[,list(chr, pos_hg38)], fromLast = F)) & (!is.na(rsid_new)), rsid_mrg:=rsid_new]
matchtable[!(duplicated(matchtable[,list(chr, pos_hg38)], fromLast = T) | duplicated(matchtable[,list(chr, pos_hg38)], fromLast = F)) & (!is.na(rsid_new)), keep:=TRUE]

# ID not listed in dbSNP (see restriction of variants in the R-package...) -> assign ID given in HapMap3 (might be old...)
matchtable[is.na(rsid_new),rsid_mrg:=rsid_old]
matchtable[is.na(rsid_new),keep:=TRUE]   

# remove rows not marked with "keep", and drop unnecessary columns
cat(paste0('Removing ',sum(!matchtable$keep),' matches because of conflicts.\n'))

setnames(matchtable, colnames(matchtable), gsub('_old$',paste0('_dbSNP',opt$dbsnp_version_hg37),colnames(matchtable)))
setnames(matchtable, colnames(matchtable), gsub('_new$',paste0('_dbSNP',opt$dbsnp_version_hg38),colnames(matchtable)))

fwrite(matchtable[!(keep)], gsub('\\.tsv(\\.gz)?$','_removed.tsv.gz', opt$out), sep='\t', scipen=50, na='.', quote=FALSE)

matchtable <- matchtable[(keep)]
matchtable[,keep:=NULL]
                               
# export
fwrite(matchtable, opt$out, sep='\t', scipen=50, na='.', quote=FALSE)