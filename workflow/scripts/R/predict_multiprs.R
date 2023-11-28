#!/usr/bin/env Rscript


start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--model_rds", action="store", default=NA, type='character',
            help="Path to saved GenoPred glmnet-model"),
make_option("--prs_dir", action="store", default="prs/", type="character")
)

opt = parse_args(OptionParser(option_list=option_list))

out_prefix <- gsub('\\model.rds$','',opt$model_rds)

if (!endsWith(opt$prs_dir, '/')){
    opt$prs_dir <- paste0(opt$prs_dir,'/')    
}

stopifnot(dir.exists(opt$prs_dir))

suppressMessages(library(data.table))
suppressMessages(library(glmnet))
suppressMessages(library(stringr))

# loading the model at the optimal hyperparameter values
m <- readRDS(opt$model_rds)
c <- coef(m, m$tuneValue$lambda)
non_zero_coef <- c[2:nrow(c),]
non_zero_coef <- non_zero_coef[non_zero_coef != 0]

cat(paste0(length(non_zero_coef), ' non-zero coefficients in model (out of ', nrow(c)-1, ')\n'))

# remove x/y to save memory
m$call$x <- NULL
m$call$y <- NULL


# identifying parameter and input files
model_rds <- gsub('^\\./','',opt$model_rds)
study_id <- strsplit(model_rds,split = '/')[[1]][4]
cat(paste0('Study ID: ',study_id,'\n'))

ancestry <- strsplit(model_rds,split = '/')[[1]][5]
cat(paste0('Ancestry: ', ancestry,'\n'))

phenotype <- strsplit(basename(model_rds),split = '\\.')[[1]][2]
cat(paste0('Phenotype: ', phenotype, '\n'))

predictor_groups_file <- paste0(dirname(model_rds),'/',study_id,'.AllMethodComp.predictor_groups')
cat(paste0('Using predictor groups file: ', predictor_groups_file, '\n'))
predictor_groups <- fread(predictor_groups_file)

read_profiles <- function(path) {
    # function that tries to "rescue" .profiles saved in the old format...
    out <- tryCatch(
        {
            header <- colnames(fread(path, sep=' ', header=T, nrow=0))
            if (all(c('FID','IID') %in% header)){
                coltypes <- list('character'=c('FID','IID'))
            } else {
                coltypes <- list('character'=c('IID'))
            }
            # this will raise a warning if fread discards the header because of an apparent "mismatch"...
            fread(path, colClasses=coltypes, sep=' ', na.strings=c('NA','.','Inf','-Inf','inf','-inf', ''), header=T, fill=FALSE, strip.white=FALSE)
        },
        warning=function(warn){
            # the file is probably space-delimited and has empty strings without quotes as NA-values:
            message(paste0("A warning was raised for file ", path,". It might contain missing values not denoted by 'NA' or '.', or not have 'IID' (and 'FID') in the header...\nFalling back to read.table (strip.white = FALSE, sep=' ')"))
            suppressWarnings(fread(text = 'n\n',nrows = 1)) # avoid warning at next call to fread
            sink(file = paste(opt$out,'.log',sep=''), append = T)
            cat('Warning: The file ',path,' has an unexpected format.\n',sep='')
            sink()
            coltypes <- c('character', 'character')
            names(coltypes) <- c('FID','IID')
            profiles <- data.table(read.table(path, colClasses=coltypes, strip.white = FALSE, sep=' ', header=T))
            if (!all(profiles[,lapply(.SD, function(x){is.numeric(x) | all(is.na(x))}),.SDcols=colnames(profiles)[-1:-2]])){
                stop(paste0('The data in ', path, " failed format checks. Make sure to use 'FID' and 'IID' as the first two columns in the header, and avoid using empty strings to indicate NA values!"))
            }
            return(profiles)
        }
    )    
    return(out)
}

predictor_groups$group <- gsub("[[:punct:]]", ".",predictor_groups$group)


pgs <- lapply( 1:nrow(predictor_groups), function(i){
        
    group <- predictor_groups$group[i]
    profiles_file <- predictor_groups$predictors[i]
    
    selected <- grep(paste0('Group_',group), names(non_zero_coef), value=TRUE)
    
    if (length(selected) == 0){
        cat('No scores selected for group',group)
        cat(', skipping...\n')
        return(NULL)
    } else {
        cat('Selected', length(selected), 'scores for group',group,'\n')    
    }
    
    profiles <- read_profiles(profiles_file)
    matched <- match(gsub('^.*PredFile[0-9]+\\.','',selected), gsub('e[\\+\\-]','e.',colnames(profiles)))
    
    if (any(is.na(matched))){
        stop('Unable to match coefficients to scoring file.') 
    }
    
    fid_iid <- profiles[,list(FID, IID)]
    pgs <- as.matrix(profiles[,..matched]) %*% non_zero_coef[selected]
    pgs <- cbind(fid_iid, pgs)
    
    coef <- data.frame(weight=non_zero_coef[selected])
    coef$variable <- colnames(profiles)[matched]
    coef$group <- group
    
    return(list(pgs=pgs, coef=coef))
    
})

coef <- rbindlist(lapply(pgs, function(x){x$coef}))

pgs <- rbindlist(lapply(pgs, function(x){x$pgs}))
pgs <- pgs[,sum(V1),by=list(FID,IID)]

fwrite(coef, paste0(out_prefix,'multiprs_weights.tsv'))
fwrite(pgs, paste0(out_prefix,'multiprs_profiles.tsv'))
system(paste0('gzip ',paste0(out_prefix,'multiprs_profiles.tsv')))