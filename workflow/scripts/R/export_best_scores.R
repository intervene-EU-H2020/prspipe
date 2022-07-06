#!/usr/bin/env Rscript

start.time <- Sys.time()

# Author: Remo Monti, 5/5/2022

# script to filter down the .score.gz and .scale files generated for each method to include only the highest performing scores
# useful when sharing data across biobanks

suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library("optparse"))

option_list = list(
make_option("--groups", action="store", type='character',
            help="Path to predictor_groups file (required)", default=NA),
make_option("--metric", action="store", default="CrossVal_R", help='Metric to use to sort models from best to worst, default: "%default"'),
make_option("--decreasing", action="store", default=TRUE, type='logical',
            help="metric ordering (if larger = better, then should be TRUE, else FALSE)"),
make_option("--N_max", action="store", type='numeric', default=20,
            help="number of best models to select for every ancestry. Total number of selected models can be greater than this."),
make_option("--ancestries", action="store", type="character", default="AFR,AMR,EAS,EUR,SAS",
            help="comma-separated list of ancestries (superpopulations) to expect, e.g., \"EUR,AMR\", default: \"%default\""),
make_option("--phenotype", default='all', type='character', 
            help='comma-separated list of phenotypes to take into account (each considered separately) or "all" to take into account all phenotypes')
)

opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$groups)){
    stop('--groups is required')
}

cat(
'### Exporting scores with highest performance ###
Analysis started at: ',as.character(start.time),'

options\n'
)
print(opt)



# setting parameters:

metric <- opt$metric
decreasing <- opt$decreasing
N_max <- max(1, round(opt$N_max, 0))
groups <- fread(opt$groups, sep=' ', header=T)
ancestries <- strsplit(opt$ancestries, ',')[[1]]
study <- strsplit(gsub('.*prs/','',groups$predictors[1]),'/',fixed = T)[[1]][2]
cat('Study: ',study, '\n')
bbid <- strsplit(gsub('.?results/','',groups$predictors[1]),'/',fixed = T)[[1]][1]
cat('Target data:       ',bbid,'\n\n')
requested_phenotypes <- strsplit(opt$phenotype, split=',')[[1]]

ancestry_pattern <- paste0('(',paste(ancestries, collapse='|'),')')

available_phenotypes <- list.files(paste0('results/',bbid,'/PRS_evaluation/',study), pattern=paste0(ancestry_pattern,'\\.AllMethodComp\\.pred_eval\\.txt'), full.names=T, recursive=T)
available_phenotypes <- unique(sapply(available_phenotypes, function(x){strsplit(basename(x),split='\\.')[[1]][2]}))

if (opt$phenotype != 'all'){
    if (!all(requested_phenotypes %in% available_phenotypes)){
        stop('Error: not all requested phenotypes are available.')
    }
    phenotypes <- requested_phenotypes
} else {
    phenotypes <- available_phenotypes
}

# helper functions:
get_megaprs_pseudoval_model <- function(path){
    
    # get the pattern to match the MegaPRS pseudovalidation model
    
    megaprs_pseudoval_path <- list.files(path, pattern = '\\.pseudo.+\\.txt$', full.names = T, recursive = F)
    
     if (length(megaprs_pseudoval_path)){
    
        cat('Guessed MegaPRS pseudovalidation log path : ',megaprs_pseudoval_path,'\n')
        cat('Reading best pseudovalidation parameters...\n')
        
        pseudoval_params<-fread(megaprs_pseudoval_path, header=F)
        
        pat <- paste0('ldak.Model', gsub('Score_','',pseudoval_params[which.max(V2)]$V1))
         
        return(pat)
         
    } else {
         
        cat('Warning: Could not guess MegaPRS log file path\n')
        cat('         Can\'t determine MegaPRS pseudovalidation model.\n')
         
        return(NA)
    }
}


get_lassosum_pseudoval_model <- function(path){
    
    # get the pattern to match the lassosum pseudovalidation model
    
    lassosum_log_path <- list.files(path, pattern = '\\.log$', full.names = T, recursive = F)
    
    if (length(lassosum_log_path) == 1){
    
        cat('Guessed lassosum log path : ',lassosum_log_path,'\n')
        cat('Reading best pseudovalidation parameters...\n')
        
        pseudoval_params<-system2(command='grep',args=c('-A','2','"Pseudovalidated parameters:"',lassosum_log_path),stdout = T, wait=T)
        
        stopifnot(startsWith(pseudoval_params[2],'s ='))
        stopifnot(startsWith(pseudoval_params[3],'lambda ='))
        
        lassosum_s <- gsub('s = ',replacement = 's',pseudoval_params[2],fixed = T)
        lassosum_lambda <- gsub('lambda = ',replacement = 'lambda',pseudoval_params[3],fixed = T)
        
        if (nchar(lassosum_lambda) > 11){
            lassosum_lambda <- gsub('[0-9]$','[0-9]+$',lassosum_lambda) 
        } else {
            lassosum_lambda <- paste0(lassosum_lambda, '$')
        }
        
        pat <- paste0(lassosum_s,'.',lassosum_lambda)
        
        return(pat)
        
    } else {
        cat('Warning: Could not guess lassosum log file path\n')
        cat('         Can\'t determine lassosum pseudovalidation model.\n')
        
        return(NA)
    }
}


# processing loop:

for (g in seq_along(groups$group)){
    
    group_name <- groups$group[g]
    
    cat('*** processing group',group_name,'***\n\n')
    
    # columns from score file to select
    colid <- numeric()
    
    # score file path 
    score_file <- paste0('prs/',group_name,'/',study,'/1KGPhase3.w_hm3.',study,'.score.gz')
    score_dir <- dirname(score_file)
    
    # Model identifiers do not contain "SCORE_" or special characters, so we get rid of them here:
    score_columns_orig <- colnames(fread(cmd=paste0('zcat ',score_file,' | head -n 1'), nrow=0, header=T))
    score_columns <- gsub('^SCORE_','',score_columns_orig)
    score_columns <- gsub('[[:punct:]]','.',score_columns)
    
    cat('processing ancestries: ')
    
    skip <- FALSE
    
    # keep track which phenotypes contributed
    phenos_processed <- list()
    
    for (ancestry in ancestries){
        
        if (skip){
            cat(paste0(ancestry,'(skipped); '))
            next
        }
        
        
        for (phenotype in phenotypes){
        
            infile <- paste0('results/',bbid,'/PRS_evaluation/',study,'/',ancestry,'/',study,'.',phenotype,'.',ancestry,'.AllMethodComp.pred_eval.txt')
            
            if (!file.exists(infile)){
                cat(paste0(ancestry,'(missing evaluation file for "',phenotype,'"); '))
                next
            }
            
            tmp_pred_eval <- fread(infile, header=T)
            
            if (nrow(tmp_pred_eval) == 0){
                cat(paste0(ancestry,'(evaluation file empty for "',phenotype,'"); '))
                next
            }
            
            # tmp_pred_eval <- tmp_pred_eval[grepl(paste0('PredFile', g, '\\.'), Model)][order(get(metric), decreasing = decreasing)]
            grp_pattern <- paste0('^',gsub('[[:punct:]]','\\.',group_name),'\\.PredFile')
            tmp_pred_eval <- tmp_pred_eval[grepl(grp_pattern, Model)][order(get(metric), decreasing = decreasing)]
            
            if (nrow(tmp_pred_eval) <= N_max){
                cat(paste0('\nNumber of scores evaluated for ', group_name,' smaller or equal to N_max=', N_max,'. Will carry over all scores.\n'))
                colid <- 1:nrow(tmp_pred_eval) + 2
                skip <- TRUE
                break
            }
            
            tmp_pred_eval <- tmp_pred_eval[1:min(nrow(tmp_pred_eval), N_max)]
            
            # the mapping from Model to score column is a mess.
            # score file colums do not contain "...PredFile..." or "_group", so we get rid of them here:
            models_select <- gsub('_group$','',gsub(paste0('.*PredFile[0-9]+\\.'),'',tmp_pred_eval$Model))
            models_select <- sapply(models_select, function(x){paste(strsplit(x,split = '.', fixed = TRUE)[[1]][-1], collapse='.')})
                    
            colid <- union(colid, match(models_select, score_columns))
            
            if (any(is.na(colid))) {
                stop(paste('Could not identify all columns for "', group_name, '" corresponding to selected Models in',infile))    
            }
                    
            cat(paste0(ancestry,'; '))
            
            if (length(phenos_processed[[ancestry]]) > 0){
                phenos_processed[[ancestry]] <- union(phenos_processed[[ancestry]], phenotype)
            } else {
                phenos_processed[[ancestry]] <- phenotype
            }
                        
        }

    }
                
    cat ('... done.\n')
    
    if (length(phenos_processed) > 0){
        cat('Phenotypes taken into account for each ancestry:\n')
        for (i in seq_along(phenos_processed)){
            p <- paste(phenos_processed[[i]], collapse=', ')
            cat(paste0(names(phenos_processed)[i],': ',p,'\n'))
        }
        cat('\n')
    }
    
    if (group_name %in% c('lassosum','ldpred2','prscs','megaprs')){
        
        # these Methods do pseudovalidation to guess a good set of hyperparameters
        # we want to always inlcude the pseudovalidation model when we share scores with others!
        
        # DBSLMM and SBayesR also do pseudovalidation but they only produce 1 model so they don't need special treatment...
        
        cat(paste0('Extracting pseudovalidation model for ', group_name, '\n'))
        
        if (group_name == 'lassosum'){
            pat <- get_lassosum_pseudoval_model(score_dir)
        } else if (group_name == 'megaprs'){
            pat <- get_megaprs_pseudoval_model(score_dir)
        } else if (group_name == 'ldpred2'){
            pat <- '.*auto'
        } else if (group_name == 'prscs'){
            pat <- '.*phi\\.auto'
        }
        
        if (!is.na(pat)){
            colid_pseudoval <- grep(pat, score_columns)
            if (length(colid_pseudoval) == 1){
                if (colid_pseudoval %in% colid){
                    cat('pseudoval model for method',group_name,' has already been selected as one of the best models.\n')
                } else {
                    cat('pseudoval model for method',group_name,' added (+1 score).\n')
                }
                colid <- union(colid, colid_pseudoval)
            } else {
                cat('Warning: Could not determine ', group_name, ' pseudoval model.\n')
            }
        }
        
    }
    
    colid_str <- paste(c(1,2,colid), collapse=',')
    
    outfile <- paste0('temp/', score_file)
    
    if (file.exists(outfile)){
        file.remove(outfile)
    }
    
    out_dir <- dirname(outfile)
    
    system(paste0('mkdir -p ', out_dir))
    system(paste0('zcat ', score_file,' |  cut -d " " -f',colid_str,' | gzip > ',outfile))
    
    if (!file.exists(outfile)){
        stop(paste0('Error extracting scores from ', score_file))
    }
    
    scale_rows_get <- score_columns_orig[colid]
    
    scale_files <- paste0(dirname(score_file),'/1KGPhase3.w_hm3.',study,'.',ancestries,'.scale')
    scalings <- lapply(scale_files, fread, header=T)
                
    for (s in seq_along(scale_files)){
        if (!(all(scale_rows_get %in% scalings[[s]]$Param))){
            stop(paste0('Missing scales in ', scale_files[s]))
        }
        scalings[[s]] <- scalings[[s]][Param %in% scale_rows_get]
        fwrite(scalings[[s]], paste0('temp/',scale_files[s]),sep=' ', col.names=T)
    }
                
    cat(paste0('Extracted best ', length(colid), ' scores from ',score_file,' to ',outfile,'\ngroup "',group_name,'" done! \n\n\n'))
    
}