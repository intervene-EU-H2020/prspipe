
library(data.table)

infiles_profiles <- snakemake@input[['pgs']]
infiles_weights <- snakemake@input[['weights']]

parse_wildcards <- function(files){
    
    splt <- strsplit(files, '/')
    
    study <- sapply(splt, '[[', 4)
    ancestry <- sapply(splt, '[[', 5)
    phenotype <- sapply(splt, function(x){strsplit(x[length(x)],'\\.')[[1]][2]})
    method <- sapply(splt, function(x){gsub('_group.*$','',gsub('^.*AllMethodComp\\.','',x[length(x)]))})
    
    metadata <- data.table(path=files, study=study, ancestry=ancestry, phenotype=phenotype, method=method)
    
    return(metadata)
}

metadata <- merge(parse_wildcards(infiles_profiles), parse_wildcards(infiles_weights), by = c('study','ancestry','phenotype','method'), suffixes = c('_p','_w'))

cat('Processing',nrow(metadata),'input files\n')

spa <- unique(metadata[,list(study, phenotype, ancestry)])

for (i in 1:nrow(spa)){

    # merge the profiles.
    
    r <- spa[i]
    
    cat(paste0('processing study=', r$study, ' phenotype=', r$phenotype, ' ancestry=', r$ancestry, '\n'))
    
    tmpmd <- metadata[study==r$study & phenotype == r$phenotype & ancestry == r$ancestry]    
    tmpdat <- lapply(tmpmd$path_p, fread)
    
    for (j in seq_along(tmpdat)){
        tmpdat[[j]]$study <- r$study
        tmpdat[[j]]$phenotype <- r$phenotype
        tmpdat[[j]]$ancestry <- paste0(r$ancestry,'.MultiPRS')
        tmpdat[[j]]$method <- tmpmd$method[j]
    }
    
    tmpdat <- dcast(rbindlist(tmpdat),FID + IID ~ study + phenotype + method + ancestry, value.var = 'V1')
        
    outfile <- paste0(dirname(tmpmd$path_p[1]), '/', r$study, '.', r$phenotype, '.', r$ancestry, '.MultiPRS_profiles.merged.txt')
    
    fwrite(tmpdat, file = outfile, sep=' ')
    system(paste0('gzip ', outfile))
    
    # merge the coef files.
    
    tmpdat <- lapply(tmpmd$path_w, fread)

    for (j in seq_along(tmpdat)){
        tmpdat[[j]]$study <- r$study
        tmpdat[[j]]$phenotype <- r$phenotype
        tmpdat[[j]]$ancestry <- r$ancestry
        tmpdat[[j]]$method <- paste0(tmpmd$method[j],'.MultiPRS')
    }
    
    coef_merged <- rbindlist(tmpdat)

    outfile <- paste0(dirname(tmpmd$path_p[1]), '/', r$study, '.', r$phenotype, '.', r$ancestry, '.MultiPRS_weights.merged.tsv')
    fwrite(coef_merged, file = outfile, sep='\t')
    system(paste0('gzip ', outfile))
    
}
