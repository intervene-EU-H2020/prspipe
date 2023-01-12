#!/usr/bin/env Rscript

library(data.table)
library(bigsnpr)
library(dplyr)

in_prefix <- commandArgs(trailingOnly=TRUE)[1]

out_prefix <- gsub('_chr$', '', in_prefix)
cat(out_prefix,'\n')

acount <- paste0(in_prefix,1:22,'.acount')
vmis <- paste0(in_prefix,1:22,'.vmiss')

afreq <- lapply(seq_along(acount), function(x){
    a <- fread(acount[x])
    v <- fread(vmis[x])
    merge(a, v, by = c('#CHROM','ID'), suffixes = c('_acount','_vmiss'))
    }
)

afreq <- rbindlist(afreq)

afreq[,ALT_freq:=ALT_CTS / OBS_CT_acount]

setnames(afreq, c("#CHROM", "ID", "REF", "ALT", "ALT_freq"), c("chr", "rsid", "a0", "a1", "freq"))
target <- afreq[,list(chr, rsid, a0, a1, freq)]

target$beta <- 1

rsid_pos <- fread('resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping_cached.tsv.gz', select=c(2,5), header = T)

target[,pos:=rsid_pos$pos_hg19[match(rsid, rsid_pos$rsid)]]

all_freq <- bigreadr::fread2('resources/bigsnpr/ref_freqs.1kg.w_hm3.tsv.gz')
projection <- bigreadr::fread2('resources/bigsnpr/projection.1kg.w_hm3.tsv.gz')

matched <- snp_match(target, all_freq[,1:5])
matched$avg_freq <- rowMeans(all_freq[matched$`_NUM_ID_`,-(1:5)])

matched <- matched[!is.na(matched$freq), ]

if (all(matched$beta == 1)){
    c_match <- cor(matched[matched$beta == 1, 'freq'], matched[matched$beta == 1, 'avg_freq'], use = 'complete')
    if (c_match < 0){
        matched$freq <- 1 - matched$freq    
    }
} else if (all(matched$beta == -1)){
    c_flipped <- cor(matched[matched$beta == -1, 'freq'], matched[matched$beta == -1, 'avg_freq'], use = 'complete')
    if (c_flipped < 0){
        matched$freq <- 1 - matched$freq    
    }
} else {
    
    c_match <- cor(matched[matched$beta == 1, 'freq'], matched[matched$beta == 1, 'avg_freq'], use = 'complete')
    c_flipped <- cor(matched[matched$beta == -1, 'freq'], matched[matched$beta == -1, 'avg_freq'], use = 'complete')
    
    if ( (c_match > 0) & (c_flipped < 0) ){
        matched$freq <- ifelse(matched$beta == 1, matched$freq, 1 - matched$freq)    
    } else if ((c_match < 0) & (c_flipped > 0)) {
        matched$freq <- ifelse(matched$beta == 1, matched$freq, 1 - matched$freq)   
    } else {
        cat('Warning: could not safely determine direction to flip alleles.\n')
        if (c_match > 0 ){
            matched$freq <- ifelse(matched$beta == 1, matched$freq, 1 - matched$freq)     
        } else {
            matched$freq <- ifelse(matched$beta == 1, matched$freq, 1 - matched$freq)   
        }
    }
}

png(paste0(out_prefix,'.freq_scatter.png'))
plot(matched$freq[1:10000], matched$avg_freq[1:10000], col=ifelse(matched$beta[1:10000]==1, 'red', 'blue'), xlab = 'matched freq', ylab = 'avg bigsnpr freq')
dev.off()

correction <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)

res <- snp_ancestry_summary(
  freq = matched$freq,
  info_freq_ref = all_freq[matched$`_NUM_ID_`, -(1:5)],
  projection = projection[matched$`_NUM_ID_`, -(1:5)],
  correction = correction
)

res <- data.frame(res)
colnames(res) <- 'prop'
res$ancestry <- row.names(res)
res <- res[c('ancestry','prop')]

write.table(res, paste0(out_prefix,'.bigsnpr_ancestry_prop_results.tsv'), sep='\t', row.names=F)

