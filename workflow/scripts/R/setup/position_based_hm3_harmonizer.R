#!/usr/bin/env Rscript

suppressMessages(library("optparse"))
option_list = list(
    make_option("--bim", action="store", default=NA, type='character', help="Path to bim file"),
    make_option("--out_prefix", action="store", default=NA, type='character', help='Output file prefix'),
    make_option("--genome_build", action="store", default="infer", type="character", help="Genome build, one of [hg19, hg38, infer], default: infer"),
    make_option("--mapping", action="store", default="resources/hapmap3/hapmap3_mapping.tsv.gz", type="character", help="mapping file."),
    make_option("--plink2", action="store", default="plink2", type="character", help="path to plink2 binary")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

if(!dir.exists(dirname(opt$out_prefix))){
    dir.create(dirname(opt$out_prefix))
}

bim = fread(opt$bim, sep = '\t', col.names = c('chr','rsid','dst','pos','a1','a2'))

# add IUPAC codes
bim[,IUPAC:='']
bim[ (a1 == 'A' & a2 =='T') | (a1 == 'T' & a2 =='A'), IUPAC:='W']
bim[ (a1 == 'C' & a2 =='G') | (a1 == 'G' & a2 =='C'), IUPAC:='S']
bim[ (a1 == 'A' & a2 =='G') | (a1 == 'G' & a2 =='A'), IUPAC:='R']
bim[ (a1 == 'C' & a2 =='T') | (a1 == 'T' & a2 =='C'), IUPAC:='Y']
bim[ (a1 == 'G' & a2 =='T') | (a1 == 'T' & a2 =='G'), IUPAC:='K']
bim[ (a1 == 'A' & a2 =='C') | (a1 == 'C' & a2 =='A'), IUPAC:='M']

mapping <- fread(opt$mapping, sep='\t')

if (opt$genome_build == 'infer'){
    l1 = nrow(merge(bim[,list(chr,pos,IUPAC)], mapping[,list(chr, pos_hg19, IUPAC)], by.x = c('chr','pos','IUPAC'), by.y=c('chr','pos_hg19','IUPAC'), suffixes=c('_target', '_ref'), all.x = FALSE, all.y = FALSE))
    l2 = nrow(merge(bim[,list(chr,pos,IUPAC)], mapping[,list(chr, pos_hg38, IUPAC)], by.x = c('chr','pos','IUPAC'), by.y=c('chr','pos_hg38','IUPAC'), suffixes=c('_target', '_ref'), all.x = FALSE, all.y = FALSE))
    if (l1 > l2){
        opt$genome_build <- 'hg19'
    } else {
        opt$genome_build <- 'hg38'
    }
    cat('Genome-build inferred to be', opt$genome_build)
}

merged = merge(bim, mapping, by.x = c('chr','pos','IUPAC'), by.y=c('chr',paste0('pos_', opt$genome_build),'IUPAC'), suffixes=c('_target', '_ref'), all.x = TRUE, no.dups = TRUE)

stopifnot(all(merged$pos == bim$pos)) # sanity check, we must keep the same order!

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

rc <- system2('bin/plink2', c('--bfile', gsub('\\.bed$','',tmp_bed),'--extract',tmp_keepfile,'--make-bed','--out',opt$out_prefix))

if (rc != 0){
    file.remove(tmp_fam, tmp_bed, tmp_bim, tmp_keepfile)
    stop('Error: call to plink2 had exit status',rc)
} else {
    file.remove(tmp_fam, tmp_bed, tmp_bim, tmp_keepfile)
}