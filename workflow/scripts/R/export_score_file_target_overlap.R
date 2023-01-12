#!/usr/bin/env Rscript

library(data.table)
library(bigsnpr)
library(dplyr)

ref <- fread('resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping_cached.tsv.gz', select=c(1:5,15), header = T)
setnames(ref, c('chr','rsid','a1','a0','pos','freq_1kg'))

read_af_biobank <- function(in_prefix){

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
        
    setnames(afreq, c("#CHROM", "ID", "REF", "ALT", "ALT_freq", "F_MISS"), c("chr", "rsid", "a0", "a1", "freq", "fmiss"))
    target <- afreq[,list(chr, rsid, a0, a1, freq, fmiss)]
    
    target$beta <- 1
    target$pos <- ref$pos[match(target$rsid, ref$rsid)]
    
    matched <- snp_match(target, ref)
    matched$freq <- with(matched, ifelse(beta == 1, freq, 1-freq))
    
    stopifnot(cor(matched$freq, matched$freq_1kg) > 0.5)
    stopifnot(all(matched$`rsid.ss` == matched$`rsid`))
    
    matched$freq_1kg <- NULL
    matched$`rsid.ss` <- NULL

    return(matched)

}


prefix <- commandArgs(trailingOnly=TRUE)[1]

out_prefix <- gsub('_chr$','',prefix)

af_miss <- data.table(read_af_biobank(prefix))

rsid_scores <- readRDS('temp/prs_rsid.rds')

refid_scores <- lapply(rsid_scores, function(x){match(x, ref$rsid)})

n_variants <- sapply(refid_scores, length)


n_overlap <- lapply(refid_scores, function(x){
   
   matched <- x %in% af_miss$`_NUM_ID_`
    
   result <- data.frame(n_overlap=sum(matched))
   result$n_missing <- sum(!matched)
    
   matched <- af_miss$`_NUM_ID_` %in% x
    
   tmp <- af_miss[(fmiss >= 0.1) & matched,list(fmiss, rsid)]

   result$n_miss_01 <- nrow(tmp)
   result$n_miss_02 <- sum(tmp$fmiss >= 0.2)
   result$n_miss_03 <- sum(tmp$fmiss >= 0.3)
   result$n_miss_04 <- sum(tmp$fmiss >= 0.4)
   result$n_miss_05 <- sum(tmp$fmiss >= 0.5)
   
    
   tmp <- af_miss[ ((freq >= 0.999) | (freq <= 0.001))  & matched, list(fmiss, freq, rsid)]

   result$n_rare <- nrow(tmp)
   result$n_invariant <- sum( (tmp$freq == 1) | (tmp$freq == 0 ) )

   return(result)
})


result <- do.call(rbind,n_overlap)

stopifnot(all(row.names(result) == names(refid_scores)))

result$n_total <- n_variants

outfile <- paste0(out_prefix,'.score_overlaps.tsv')

write.table(result, file=outfile, col.names = NA, sep='\t')

system(paste0('gzip ', outfile))



