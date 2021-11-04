#!/usr/bin/env Rscript

suppressMessages(library("optparse"))
option_list = list(
    make_option("--bim", action="store", default=NA, type='character', help="Path to bim file from 1000 Genomes"),
    make_option("--hm3_bim", action="store", default=NA, type="character", help="Path to bim file from HapMap3"),
    make_option("--chr", action="store", default=NA, type="character", help="chromosome"),
    make_option("--dbsnp_version", action="store", default="144", type="character"),
    make_option("--out", action="store", default=NA, type='character', help='Output file.'),
    make_option("--plink2", action="store", default="plink2", type="character", help="path to plink2 binary")
)

opt = parse_args(OptionParser(option_list=option_list))

suppressMessages(library(data.table))
# hg19 bim
bim = opt$bim

# the dbSNP version
dbsnp_version = opt$dbsnp_version
# chromosome
chrom = opt$chr
# hapmap3 bim(s)
hm3 = strsplit(opt$hm3_bim,split=',')[[1]]

# load dbSNP from the custom location
.libPaths(c(.libPaths(), "R/local/"))
suppressMessages(library('GenomicRanges'))
suppressMessages(library(paste0("SNPlocs.Hsapiens.dbSNP",dbsnp_version,".GRCh37"),character.only = TRUE))

bim <- fread(bim, sep = '\t', col.names = c('chr','rsid','dst','pos','a1','a2'))

# add IUPAC codes
bim[,IUPAC:='.']
bim[ (a1 == 'A' & a2 =='T') | (a1 == 'T' & a2 =='A'), IUPAC:='W']
bim[ (a1 == 'C' & a2 =='G') | (a1 == 'G' & a2 =='C'), IUPAC:='S']
bim[ (a1 == 'A' & a2 =='G') | (a1 == 'G' & a2 =='A'), IUPAC:='R']
bim[ (a1 == 'C' & a2 =='T') | (a1 == 'T' & a2 =='C'), IUPAC:='Y']
bim[ (a1 == 'G' & a2 =='T') | (a1 == 'T' & a2 =='G'), IUPAC:='K']
bim[ (a1 == 'A' & a2 =='C') | (a1 == 'C' & a2 =='A'), IUPAC:='M']

options(scipen = 999)

snp_pos <- GRanges(seqnames=bim$chr, ranges=IRanges(start=bim$pos, width=1), rsid=bim$rsid)

cat(paste0('Overlapping ',length(snp_pos),' SNP positions with dbSNP annotation ',dbsnp_version,'.\n'))

# get the overlaps by position with dbSNP
snp_overlaps <- snpsByOverlaps(get(paste0('SNPlocs.Hsapiens.dbSNP',dbsnp_version,'.GRCh37')), ranges = snp_pos)
cat(paste0('Found ', length(snp_overlaps),' (',round(length(snp_overlaps)/length(snp_pos)*100,2), '%) SNPs with rsIDs in dbSNP.\n'))

snp_overlaps<-data.table(as.data.frame(snp_overlaps))
setnames(snp_overlaps, old=c('seqnames','pos','strand','RefSNP_id','alleles_as_ambig'), new=c('chr','pos','strand','rsid','IUPAC'))
snp_overlaps[,strand:=NULL]
bim$chr <- as.character(bim$chr)
snp_overlaps$chr <- as.character(snp_overlaps$chr)

bim[,i_bim:=seq(1,nrow(bim))]

bim_mrg <- merge(bim, snp_overlaps, by = c('chr','pos'), suffixes = c('_1kg','_dbSNP'), all.x = TRUE)
stopifnot(nrow(bim_mrg) > 0) # if the intersection is empty. Something is wrong...

hm3_snps <- lapply(hm3, function(x){fread(cmd=paste0("awk '$1 ==  ",chrom,"' ",x), select = 2, header=F)})
hm3_snps <- unique(rbindlist(hm3_snps))

cat(paste0('Overlapping with ',nrow(hm3_snps), ' HapMap3 variants\n'))

bim_mrg[,in_hm3 := (rsid_dbSNP %in% hm3_snps$V2)]
bim_mrg[,sum_in_hm3:=sum(in_hm3),by=list(i_bim)]
bim_mrg[,is_duplicated:=FALSE]
bim_mrg[duplicated(i_bim)|duplicated(i_bim,fromLast=T), is_duplicated:=TRUE]
bim_mrg[,drop:=FALSE]
bim_mrg[(sum_in_hm3 == 1) & (is_duplicated) & (in_hm3 == FALSE), drop:=TRUE]

# handling of retired rs numbers with both are in hm3 -> if rs numbers are merged, the higher number is retired!
bim_mrg[is_duplicated & sum_in_hm3 > 1, rsnumber:=as.numeric(gsub('rs','',rsid_dbSNP)) ]
# drop all except the one with the lowest rs-number:
bim_mrg[is_duplicated & sum_in_hm3 > 1,drop:=rsnumber>min(rsnumber),by = list(chr, pos, rsid_1kg)]

bim_mrg = bim_mrg[!(drop)]
# does this actually do anything anymore?
bim_mrg = bim_mrg[ ((!duplicated(i_bim)) & (sum_in_hm3 == 0)) | (sum_in_hm3 >= 1) ][order(i_bim)]

stopifnot(all(bim_mrg$i_bim == bim$i_bim)) # sanity check
cat(paste0('recovered ', sum(bim_mrg$in_hm3), ' (', round(100*sum(bim_mrg$in_hm3)/nrow(hm3_snps),2),'%) of HapMap3 variants.'))

# in_fam <- gsub('\\.bim$','.fam', opt$bim)
# in_bed <- gsub('\\.bim$','.bed', opt$bim)
# 
# tmp_fam <- paste0(opt$out_prefix, '_tmp.fam')
# tmp_bed <- paste0(opt$out_prefix, '_tmp.bed')
# tmp_bim <- paste0(opt$out_prefix, '_tmp.bim')
# 
# fwrite(bim_mrg[,list(chr, rsid_dbSNP, dst, pos, a1, a2)], tmp_bim, scipen = 50, sep='\t', col.names = F, row.names = F, na = '.', quote=F)
# 
# system(paste('ln -s -r',in_fam,tmp_fam))
# system(paste('ln -s -r',in_bed,tmp_bed))
# 
# tmp_keepfile <- paste0(opt$out_prefix,'.keep')
# fwrite(bim_mrg[(in_hm3),list(rsid_dbSNP),], scipen = 50, row.names = F, col.names = F, file = tmp_keepfile)
# 
# rc <- system2(opt$plink2, c('--bfile', gsub('\\.bed$','',tmp_bed),'--extract',tmp_keepfile,'--make-bed','--out',opt$out_prefix))
# 
# if (rc != 0){
#     rmcode = file.remove(tmp_fam, tmp_bed, tmp_bim, tmp_keepfile)
#     stop('Error: call to plink2 had exit status',rc)
# } else {
#     rmcode = file.remove(tmp_fam, tmp_bed, tmp_bim, tmp_keepfile)
# }

setnames(bim_mrg, 'pos', 'pos_hg19')
fwrite(bim_mrg[(in_hm3) | ((rsid_1kg != '.') & (!is.na(rsid_1kg)) & (rsid_1kg %in% hm3_snps$V2)),list(chr, pos_hg19, rsid_1kg, rsid_dbSNP, a1, a2, IUPAC_1kg, IUPAC_dbSNP)], opt$out, sep='\t', quote=F, na=".", row.names=F)
