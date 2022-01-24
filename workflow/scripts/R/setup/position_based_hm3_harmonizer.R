#!/usr/bin/env Rscript

suppressMessages(library("optparse"))
option_list = list(
    make_option("--bim", action="store", default=NA, type='character', help="Path to bim file"),
    make_option("--out_prefix", action="store", default=NA, type='character', help='Output file prefix'),
    make_option("--genome_build", action="store", default="infer", type="character", help="Genome build, one of [hg19, hg38, infer], default: infer"),
    make_option("--mapping", action="store", default=NA, type="character", help="mapping file."),
    make_option("--plink2", action="store", default="plink2", type="character", help="path to plink2 binary"),
    make_option("--rsid_col", action="store", default="rsid_mrg", type="character"),
    make_option("--plink2_args", action="store", default="--make-bed", type="character", help="arguments passed to plink2 (use this only if you know what you are doing)")
)

opt = parse_args(OptionParser(option_list=option_list))

suppressMessages(library(data.table))

if(!dir.exists(dirname(opt$out_prefix))){
    dir.create(dirname(opt$out_prefix))
}

bim = fread(opt$bim, sep = '\t', col.names = c('chr','rsid','dst','pos','a1','a2'))

# add IUPAC codes
bim[,IUPAC:='.']
bim[ (a1 == 'A' & a2 =='T') | (a1 == 'T' & a2 =='A'), IUPAC:='W']
bim[ (a1 == 'C' & a2 =='G') | (a1 == 'G' & a2 =='C'), IUPAC:='S']
bim[ (a1 == 'A' & a2 =='G') | (a1 == 'G' & a2 =='A'), IUPAC:='R']
bim[ (a1 == 'C' & a2 =='T') | (a1 == 'T' & a2 =='C'), IUPAC:='Y']
bim[ (a1 == 'G' & a2 =='T') | (a1 == 'T' & a2 =='G'), IUPAC:='K']
bim[ (a1 == 'A' & a2 =='C') | (a1 == 'C' & a2 =='A'), IUPAC:='M']

mapping <- fread(opt$mapping, sep='\t')
nr <- nrow(mapping)
cat('Loaded ', nr, ' variants from mapping file.\n')
mapping <- mapping[ a1 %in% c('A','C','T','G') & a2 %in% c('A', 'C', 'T', 'G')  ]
cat('dropped ', nr-nrow(mapping), ' indels.\n')


if (opt$rsid_col %in% colnames(mapping)){
    setnames(mapping,opt$rsid_col,'rsid')
} else {
    stop(paste0('Error: column "',opt$rsid_col, '" not found in mapping file. available columns: ',paste(colnames(mapping), collapse=',')))
}

setnames(mapping, 'IUPAC_1kg', 'IUPAC')

if (opt$genome_build == 'infer'){
    l1 = nrow(merge(bim[,list(chr,pos,IUPAC)], mapping[,list(chr, pos_hg19, IUPAC)], by.x = c('chr','pos','IUPAC'), by.y=c('chr','pos_hg19','IUPAC'), suffixes=c('_target', '_ref'), all.x = FALSE, all.y = FALSE))
    l2 = nrow(merge(bim[,list(chr,pos,IUPAC)], mapping[,list(chr, pos_hg38, IUPAC)], by.x = c('chr','pos','IUPAC'), by.y=c('chr','pos_hg38','IUPAC'), suffixes=c('_target', '_ref'), all.x = FALSE, all.y = FALSE))
    if (l1 > l2){
        opt$genome_build <- 'hg19'
    } else {
        opt$genome_build <- 'hg38'
    }
    cat('Genome-build inferred: ', opt$genome_build, '\n')
}

bim[,i_bim:=1:nrow(bim)]

merged = merge(bim, mapping, by.x = c('chr','pos','IUPAC'), by.y=c('chr',paste0('pos_', opt$genome_build),'IUPAC'), suffixes=c('_target', '_ref'), all.x = TRUE, no.dups = TRUE)

merged <- merged[order(i_bim)]

stopifnot(all(merged$i_bim == bim$i_bim)) # sanity check, we must keep the same order!

merged[,rsid_target:=rsid_ref]
merged[is.na(rsid_target),rsid_target:='.']

in_fam <- gsub('\\.bim$','.fam', opt$bim)
in_bed <- gsub('\\.bim$','.bed', opt$bim)

tmp_fam <- paste0(opt$out_prefix, '_tmp.fam')
tmp_bed <- paste0(opt$out_prefix, '_tmp.bed')
tmp_bim <- paste0(opt$out_prefix, '_tmp.bim')

fwrite(merged[,list(chr, rsid_target, dst, pos, a1_target, a2_target)], tmp_bim, scipen = 50, sep='\t', col.names = F, row.names = F)

system(paste('ln -s -r',in_fam,tmp_fam))
system(paste('ln -s -r',in_bed,tmp_bed))

tmp_keepfile <- paste0(opt$out_prefix,'.keep')

fwrite(merged[rsid_target!='.',list(rsid_target),], scipen = 50, row.names = F, col.names = F, file = tmp_keepfile)


plink2_args <- c(c('--bfile', gsub('\\.bed$','',tmp_bed),'--extract',tmp_keepfile,'--out',opt$out_prefix), strsplit(opt$plink2_args, ' ')[[1]])
cat('Calling plink 2 with arguments: ', plink2_args, '\n')

rc <- system2(opt$plink2, plink2_args)

if (rc != 0){
    file.remove(tmp_fam, tmp_bed, tmp_bim, tmp_keepfile)
    stop('Error: call to plink2 had exit status ',rc)
} else {
    file.remove(tmp_fam, tmp_bed, tmp_bim, tmp_keepfile)
}
