#!/usr/bin/env Rscript

suppressMessages(library("optparse"))
option_list = list(
    make_option("--frq_path_prefix", action="store", default=NA, type='character', help="Path to directory containing the ancestry MAF"),
    make_option("--mapping_prefix", action="store", default=NA, type='character', help="Mapping file name prefix"),
    make_option("--rsid_col", action="store", default="rsid_mrg", type="character"),
    make_option("--min_maf", action="store", default="0.01", type="character"),
    make_option("--out_prefix", action="store", type="character")
)

# script that merges the mapping files with MAFs and filters SNPs

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

opt$frq_path_prefix <- gsub('//','/',opt$frq_path_prefix)
opt$min_maf <- as.numeric(opt$min_maf)

infiles <- list.files(opt$frq_path_prefix, pattern='chr[0-9]+\\.frq$', full.names = T, recursive = T)

mapping_files <- list.files(dirname(opt$mapping_prefix), pattern=paste0(basename(opt$mapping_prefix),'[0-9]+.tsv.gz'), full.names = T)

pop <- sapply(infiles, function(x){rev(strsplit(x, '/')[[1]])[2]})


mafs <- lapply(infiles, fread, select=c('CHR','SNP','A1','A2','MAF'))

for (i in seq_along(infiles)){
    mafs[[i]]$superpop <- paste0('MAF_',pop[i])
}

mafs <- rbindlist(mafs)

mafs <- dcast(mafs, CHR + SNP + A1 + A2 ~ superpop, value.var = c('MAF'))

mapping <- lapply(mapping_files, fread)

mapping <- rbindlist(mapping, use.names = T, fill = T)

setnames(mafs, c('CHR','SNP','A1','A2'), c('chr','rsid','a1','a2'))

setnames(mapping, opt$rsid_col, 'rsid')

mapping <- merge(mapping, mafs, by = c('chr','rsid','a1','a2'))

mapping <- mapping[order(chr, pos_hg19, a1, a2)]

setcolorder(mapping, c('chr','rsid','a1','a2','pos_hg19','pos_hg38',grep('rsid_', colnames(mapping), value = T), grep('IUPAC', colnames(mapping), value=T), grep('MAF', colnames(mapping), value=T)))

mafcols <- grep('MAF_', colnames(mapping), value=T)
mafcols <- mafcols[ mafcols != 'MAF_AllAncestry']

maf_mask <- rowSums(mapping[,mafcols,with=F][,lapply(.SD, function(x){ifelse(x<0.5, x, 1-x)})] > opt$min_maf) > 0
dup_mask <- !duplicated(mapping$rsid, fromLast=T) | duplicated(mapping$rsid, fromLast = F)
# Retain non-ambiguous
abig_mask <-!(mapping$a1 == 'A' & mapping$a2 == 'T' | 
              mapping$a1 == 'T' & mapping$a2 == 'A' | 
              mapping$a1 == 'C' & mapping$a2 == 'G' |
              mapping$a1 == 'G' & mapping$a2 == 'C')

cat('Found ', sum(!maf_mask), ' SNPs with MAF < ', opt$min_maf,' in all superpopulations. Removing...\n')
cat('Found ', sum(!dup_mask), ' SNPs with duplicated identifiers. Removing...\n')
cat('Found ', sum(!abig_mask), ' SNPs with ambiguous strand. Removing...\n')

keep <- (maf_mask & dup_mask & abig_mask)

cat('Removing ', sum(!keep), ' SNPs based on the filters above.(',paste0(opt$out_prefix, '_rm.tsv.gz'),')\n')
cat('Writing remaining ', sum(keep), ' SNPs to ', paste0(opt$out_prefix, '.tsv.gz\n'))

fwrite(mapping[!keep], paste0(opt$out_prefix, '_rm.tsv.gz'), sep='\t')
fwrite(mapping[keep], paste0(opt$out_prefix, '.tsv.gz'), sep='\t')

