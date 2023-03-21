#!/usr/bin/env Rscript


start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--target_chr", action="store", default=NA, type='character',
            help="Prefix to per-chromosome target genotype data"),
make_option("--keep", action="store", default=NA, type="character",
            help="Path to optional keep-file to subset genotype data"),
make_option("--extract", action="store", default=NA, type="character", 
            help="Path to optional file to subset variants"),
make_option("--score", action="store", type="character", 
            help="Path tp scoring file"),
make_option("--afreq", action="store", default=NA, type="character",
            help="Prefix for per-chromosome plink2 afreq files."),
make_option("--out", action="store", default=NA, type="character",
            help="Output prefix"),
make_option("--keep_pattern", action="store", default=NA, type="character",
            help="Output prefix"),
make_option("--mem", action="store", default=8000, type="numeric",
            help="Memory in Mb"),
make_option("--plink2", action="store", default="plink2", type="character",
            help="path to plink2 binary."),
make_option("--hla_only", action="store", default=FALSE, type="logical")
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)

if(!dir.exists(dirname(opt$out))){
    dir.create(dirname(opt$out), recursive = T)
}

if (endsWith(opt$score, 'gz')){
    header_cols <- strsplit(system(paste0('zcat ', opt$score, ' | head -n 1'), intern = T), split = '\t')[[1]]
} else {
    header_cols <- strsplit(system(paste0('head -n 1 ', opt$score), intern = T), split = '\t')[[1]]
}
if (length(header_cols) == 1){
    header_cols <- strsplit(header_cols, split=' ')[[1]]
    if (length(header_cols) == 1){
        stop('scoring files must be space or tab separated.')
    }
}

cat(paste0('Score file contains ', length(header_cols)-2, ' scores\n'))

ncol <- length(header_cols)

if (!is.na(opt$keep_pattern)){
    keep_columns <- grep(opt$keep_pattern, header_cols)
    if (length(keep_columns) == 0){
        stop(paste0('None of the columns in the scoring file match --keep_pattern (', opt$keep_pattern,')'))
    }
    cat(paste0('Subsetting to ', length(keep_columns), ' scores based on --keep-pattern (', opt$keep_pattern,')\n'))
    keep_columns <- paste(keep_columns, collapse = ',')
} else {
    keep_columns <- paste0('3-', ncol)
}
if (keep_columns == '3-3'){
    keep_columns <- '3'
}


profiles <- NULL

for (chrom in 1:22){
    
    if (opt$hla_only){
        if (chrom != 6){
            next
        }
        out <- paste0(opt$out, '.profiles_hla.chr', chrom)
    } else {
        out <- paste0(opt$out, '.profiles.chr', chrom)
    }
    
    cat('Scoring for chromosome', chrom, '\n')

    geno_file <- paste0(opt$target_chr, chrom)
    freq_file <- paste0(opt$afreq, chrom, '.afreq')
    
    cl_args <- c('--memory', opt$mem*0.9, '--threads', 1,'--bfile', geno_file, '--score', opt$score, 'header-read', 'cols=fid,scoresums', '--score-col-nums', keep_columns, '--out', out)
    if (!is.na(opt$keep)){
        cl_args <- c(cl_args, c('--keep', opt$keep))    
    }
    if (!is.na(opt$afreq)){
        cl_args <- c(cl_args, c('--read-freq', freq_file))    
    }
    if (!is.na(opt$extract)){
        cl_args <- c(cl_args, c('--extract', strsplit(opt$extract, split=' ')[[1]]))  
    }
    if (opt$hla_only){
        if (!is.na(opt$extract)){
            stop('Currently it is not supported to supply both --extract and --hla_only')
        } else {
            hla <- data.frame(chr=6, start=28e6, end=34e6)
            options(scipen = 999)
            write.table(hla, paste0(opt$out,'.hla.bed'), sep='\t', quote=F, col.names=F, row.names=F)
            options(scipen = 0)
            cl_args <- c(cl_args, c('--extract','bed0',paste0(opt$out,'.hla.bed')))
        }
    }
    
    # cat('plink2',paste(cl_args, collase=' '), '\n')
    
    called <- system2(opt$plink2, args = cl_args, stdout = T, stderr = T)
    
    sscore <- paste0(out, '.sscore')
    
    if (!file.exists(sscore)){
        cat(paste0('Warning: no scores produced for chromosome ',chrom,'. Check plink log file for errors.\n'))
        if (opt$hla_only){
            e <- system(paste0('grep Error: ', out, '.log'), intern=TRUE)
            if (e == 'Error: No valid variants in --score file.'){
                cat(paste0('Warning: Running with --hla_only. It is possible that no variants in the score file are in the HLA-region. Will create an empty profiles file.\n'))
                system(paste('touch ',paste0(opt$out,'.hla_profiles.tsv.gz')))
                cat('created empty file ',paste0(opt$out,'.hla_profiles.tsv.gz\n'))
                quit()
            } else {
                print(paste0('plink error ->',e))
                stop('encountered an unexpected plink error. Check plink log file for details.')
            }
        }
        next
    }
    
    if (is.null(profiles)){
        
        profiles <- fread(sscore)
        
        fid_iid <- profiles[,c(1,2)]
        profiles[,c('#FID','IID'):=NULL]
        
    } else {
        
        tmp <- fread(sscore)
        
        if (!all(tmp$IID == fid_iid$IID)){
            stop('IID mismatch while scoring.')    
        }
        
        if (!all(colnames(profiles) %in% colnames(tmp))){
            stop('score column mismatch while scoring')    
        }
        
        if (!all(colnames(tmp)[3:ncol(tmp)] %in% colnames(profiles))){
            stop('score column mismatch while scoring')    
        }
        
        tmp[,c('#FID','IID'):=NULL]
        
        for (col in colnames(profiles)) profiles[,(col):=get(col)+tmp[,col,with=FALSE]]
        
    }
    
    system(paste0('rm ', sscore))
}

profiles[,FID:=fid_iid$`#FID`]
profiles[,IID:=fid_iid$IID]

setcolorder(profiles, c(ncol(profiles)-1, ncol(profiles), 1:(ncol(profiles)-2)))

if (!opt$hla_only){
    outfile <- paste0(opt$out,'.profiles.tsv.gz')
} else {
    outfile <- paste0(opt$out,'.hla_profiles.tsv.gz')    
}

fwrite(profiles, file = outfile, sep="\t")
cat(paste0('Wrote predictions to ',outfile,'\n'))

