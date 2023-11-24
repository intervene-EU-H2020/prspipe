#!/usr/bin/env Rscript

########################################################################################
# This script can be used to evaluate existing PGS stored in different .profiles files 
# It provides more in-depth metrics compared to those produced by Model_builder_V2.R including confidence intervals  
# This script was used to generate the metrics reported in Monti et al 2023            
#########################################################################################

suppressMessages(library("optparse"))
start.time <- Sys.time()

option_list = list(
make_option("--profiles", action="store", default=NA, type='character',
            help="Path to .profiles file(s) containing the pre-calculated polygenic scores (comma-separated list)."),
make_option("--profiles_2", action="store", default=NA, type="character",
           help="Path to ukbb MultiPRS .profiles file(s) containing the pre-calculated polygenic scores (comma-separated list), columns are filtered by --keep"),
make_option("--keep", action="store", default="", type="character",
           help="Pattern used to filter the columns in profiles_2"),
make_option("--pheno", action="store", default=NA, type='character',
            help="path to phenotype file (used for subsetting)."),
make_option("--train_test", action="store", default=NA, type='character',
            help="prefix to train-test split files. (ending in test_ind.txt.gz and train_ind.txt.gz)"),
make_option("--sex", action="store", default=NA, type='character',
           help='path to file containing sex coded as 0 or 1, no header, three columns: FID, IID, SEX'),
make_option("--out", action="store", default=NA, type='character',
           help='output prefix'),
make_option("--calc_cor", action="store",default=T, type='logical',
           help='whether to calculate score-score correlations'),
make_option("--debug", action="store", default=F, type='logical',
           help='downsample number of individuals to 1000 for debugging.'),
make_option("--boot", action="store", default=0, type='numeric',
           help='number of bootstrap samples to to calculate confidence intervals for binary phenotypes R2_OBS'),
make_option("--sample_cor", action="store", default=-1, type='numeric', 
           help='unmber of individuals to sample to compute score-score correlation. default = -1 (all)'),
make_option("--threads", action="store", default=1, type='numeric',
           help='set number of cores for parallel bootstrapping')
)

suppressMessages(library(data.table))
suppressMessages(library(pROC))
suppressMessages(library(caret))
suppressMessages(library(boot))
suppressMessages(library(doMC))

set.seed(1)

printlog <- function(..., time=F){
    t <- Sys.time()
    if (!time){
        cat(paste0(paste0(...,sep=''),'\n'))
    } else {
        cat(paste0('(',as.character(t),') ',paste0(...,sep=''),'\n'))
    }
    return(invisible(t))
}


time_taken <- function(start_time=start.time){
    t <- Sys.time() - start_time
    return(paste(as.character(round(t,2)),attr(t, 'units')))
}

printlog('Analysis started at ',as.character(start.time))

opt <- parse_args(OptionParser(option_list=option_list))

printlog('Options are: ')
print(opt)

doMC::registerDoMC(cores=opt$threads)

for (i in seq_along(opt)){
    if(!is.na(opt[[i]])){
        if(opt[[i]] == 'NA') opt[[i]] <- NA
    }
}

if (!is.na(opt$train_test)){

    printlog('Reading train test splits.\n')

    train_ind <- fread(paste0(opt$train_test,'.train_ind.txt.gz'), header=F)$V1
    test_ind <- fread(paste0(opt$train_test,'.test_ind.txt.gz'), header=F)$V1
    stopifnot(length(intersect(train_ind, test_ind))==0)
}

if (!is.na(opt$sex)){

    printlog('Reading sex information.\n')

    sex <- fread(opt$sex, header=F)
    sex$fid_iid <- paste(sex$V1, sex$V2, sep=':')
    sex[,c('V1','V2'):=NULL]
    setcolorder(sex, c(2,1))
    setnames(sex, old='V3', new='sex')
}

parse_method <- function(filename){
    cat(filename,'\n')
    if (endsWith(filename, '.MultiPRS_profiles.merged.txt.gz')){
        return('MultiPRS')
    } else {
        return(gsub('[[:punct:]]','.',strsplit(filename, split='/')[[1]][4]))
    }
}

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

profile_paths <- strsplit(opt$profiles,',')[[1]]

t0 <- printlog('Reading the profile files specified by --profiles', time=T)

# reading the data

pheno <- fread(opt$pheno, header=F)
pheno <- pheno[complete.cases(pheno)]
pheno$fid_iid <-  paste0(pheno$V1, ':', pheno$V2)
pheno[,c('V1','V2'):=NULL]
setcolorder(pheno, c(2,1))
setnames(pheno, 'V3', 'Outcome')

profiles <- read_profiles(profile_paths[1])

profiles$fid_iid <- paste0(profiles$FID, ':', profiles$IID)
profiles[,c('FID','IID'):=NULL]
setcolorder(profiles, c('fid_iid', colnames(profiles)[-ncol(profiles)]))

method <- parse_method(profile_paths[1])
if (method != 'MultiPRS'){
    setnames(profiles, colnames(profiles)[-1], paste(colnames(profiles)[-1], method, sep='_'))
}

profiles <- profiles[fid_iid %in% pheno$fid_iid]

for (i in 2:length(profile_paths)){
    
    profiles_tmp <- read_profiles(profile_paths[i])
    profiles_tmp$fid_iid <- paste0(profiles_tmp$FID, ':', profiles_tmp$IID)
    profiles_tmp[,c('FID','IID'):=NULL]
    setcolorder(profiles_tmp, c('fid_iid', colnames(profiles_tmp)[-ncol(profiles_tmp)]))
    
    method <- parse_method(profile_paths[i])
    if (method != 'MultiPRS'){
        setnames(profiles_tmp, colnames(profiles_tmp)[-1], paste(colnames(profiles_tmp)[-1], method, sep='_'))
    }
    
    profiles_tmp <- profiles_tmp[fid_iid %in% pheno$fid_iid]
    profiles <- merge(profiles, profiles_tmp, by = 'fid_iid', all.x = T, all.y = T)
    
    rm(profiles_tmp)
}

printlog('finished reading data from profile files specified by --profiles, reading took ', time_taken(t0),'\n')

strip_names <- function(c){
    gsub('_[0-9a-zA-Z\\.]+_[0-9a-zA-Z\\.-]+$','',c)
}


if (!is.na(opt$profiles_2)){
    profile_paths <- strsplit(opt$profiles_2,',')[[1]]

    t0 <- printlog('Adding scores from ', length(profile_paths), ' profile files specified with --profiles_2', time=T)

    for(i in seq_along(profile_paths)){

        cat(profile_paths[i],'\n')        
        profiles_tmp <- fread(profile_paths[i], header=T)

        if (opt$keep != ""){
            stripped <- strip_names(colnames(profiles_tmp))
            matched <- which(stripped == opt$keep)
            if (length(matched) == 0){
                stop('None of the columns in profiles_2 file match pattern.')
            } else {
                cat('Found ',length(matched),'scores matching --keep pattern.\n')
            }
            matched <- c(1,2,matched)
            profiles_tmp <- profiles_tmp[,..matched]
        }
        
        profiles_tmp$fid_iid <- paste0(profiles_tmp$FID, ':', profiles_tmp$IID)
        profiles_tmp[,c('FID','IID'):=NULL]

        profiles_tmp <- profiles_tmp[fid_iid %in% pheno$fid_iid]
        profiles <- merge(profiles, profiles_tmp, by = 'fid_iid', all.x = T, all.y = T)

        rm(profiles_tmp)
    }

    printlog('Finished reading profile files specified with --profiles_2, reading took ', time_taken(t0),'\n')

}

if (all(pheno$Outcome %in% c(0,1))){
    is_binary <- TRUE
    family <- 'binomial'
    cat('Phenotype is binary.\n')
} else {
    is_binary <- FALSE
    family <- 'gaussian' # not used.
    cat('Phenotype is continuous.\n')
}

if (!is.na(opt$train_test)){
    stopifnot(all(train_ind %in% profiles$fid_iid))
    stopifnot(all(test_ind %in% profiles$fid_iid))
    profiles$split <- ''
    profiles[fid_iid %in% train_ind, split:='train']
    profiles[fid_iid %in% test_ind, split:='test']
    profiles <- profiles[!(split == '')]
} else {
    profiles[,split:='train']
}

printlog('Profiles contain ',nrow(profiles),' individuals after merging.')

if (!is.na(opt$sex)){
    profiles[match(fid_iid, sex$fid_iid),sex:=sex$sex[match(fid_iid, sex$fid_iid)]]
    rm(sex)
    if (any(is.na(profiles$sex))){
        cat('Warning: missing sex for',sum(is.na(profiles$sex)),'individuals.\nWill keep only complete cases.\n')
        profiles <- profiles[!is.na(sex)]
    }
    stopifnot(nrow(profiles)>0)
    stopifnot(all(profiles$sex %in% c(0,1)))
    non_predictor_cols <- c('fid_iid','split','sex','Outcome')
} else {
    non_predictor_cols <- c('fid_iid','split','Outcome')
}

profiles$Outcome <- pheno$Outcome[match(profiles$fid_iid, pheno$fid_iid)]
predictor_cols <- colnames(profiles)[!(colnames(profiles) %in% non_predictor_cols)]
printlog('Calculating metrics for ',length(predictor_cols),' predictors.')

if (opt$debug){
    if (1000 / nrow(profiles) > 1.){
        cat('--debug flag has no effect because number of samples is smaller than 1000.\n')
    } else {
        cat('--debug flag set, downsampling to 1000 samples.\n')
        train_index <- createDataPartition(profiles$Outcome, p = 1000/nrow(profiles), times = 1, list = F)
        profiles <- profiles[train_index]
    }
}

means <- profiles[,lapply(.SD, mean, na.rm=T),by=split,.SDcols=c(predictor_cols,'Outcome')][,stat:='Mean']
sds <- profiles[,lapply(.SD, sd, na.rm=T),by=split,.SDcols=c(predictor_cols,'Outcome')][,stat:='SD']

stats <- rbindlist(list(means, sds), use.names = T, fill = TRUE)
stats <- dcast(melt(stats, id.vars = c('split','stat')),formula = variable ~ stat + split)

rm(means, sds)
printlog('Writing means and standard deviations to ',paste0(opt$out,'.mean_sd.tsv'),'\n')
fwrite(stats, paste0(opt$out,'.mean_sd.tsv'), sep='\t', na = 'NA', quote=F)

sd_cols <- grep('SD', colnames(stats), value = T)
exclusions <- as.character(stats$variable[ rowSums(is.na(stats[,-1])) > 0 ]) # remove if any column is NA 
exclusions <- unique(c(exclusions, as.character(stats$variable[ !(stats$variable %in% exclusions) & rowSums(stats[,..sd_cols] == 0.) > 0 ]))) # remove if the standard deviation is 0.

if('Outcome' %in% exclusions){
    stop('The Outcome variable is constant or missing\n')
}

if (length(exclusions) > 0){
    cat('Warning: some predictors are constant or missing, excluding: \n', paste(exclusions),'\n\n')
    profiles[,(exclusions):=NULL]
    predictor_cols <- predictor_cols[!(predictor_cols %in% exclusions)]
}

setcolorder(profiles, c(non_predictor_cols,predictor_cols))

M <- sapply(predictor_cols, function(x){stats$Mean_train[stats$variable == x]})
SD <- sapply(predictor_cols, function(x){stats$SD_train[stats$variable == x]})

t0 <- printlog('Scaling profiles', time=T)
for (i in seq_along(predictor_cols)){
    # scale profiles by training set mean and standard deviation
    profiles[,(predictor_cols[i]):=(get(predictor_cols[i])-M[i])/SD[i]]
}
printlog('Scaling profiles took ',time_taken(t0),'\n')
rm(M,SD,sd_cols)

if (!is_binary){
    # scale outcome by training set mean and standard deviation
    profiles[,Outcome:=(Outcome - stats$Mean_train[stats$variable == 'Outcome'])/stats$SD_train[stats$variable == 'Outcome']]
}

matrix_merge <- function(x, y){
    if (is.matrix(x)){
        x <- data.table(x, keep.rownames = T)
    }
    if (is.matrix(y)){
        y <- data.table(y, keep.rownames = T)
    }
    return(data.table::merge.data.table(x,y,by = 'rn'))
}

slim_glm <- function(df){
    m <- glm(Outcome ~ ., data = df, family = family, x=F, y=F, model=F)
    m$data <- NULL
    return(m)
}

slim_lm <- function(df){
    m <- lm(Outcome ~ ., data = df, x=F, y=F, model=F)
    return(m)
}

calculate_metrics <- function(df, rescale_r2=F){
    
    # rescale_r2 should be set to TRUE if the variance of the predictors (sepcified by predictor_cols) in df is not 1.
    # predictors are not rescaled inside this function -> make sure the scales are appropriately set before (e.g. by standardization).
    
    if (rescale_r2){
        vars <- unlist(df[,lapply(.SD, var),.SDcols=predictor_cols])
        vars <- vars[predictor_cols]
    } else {
        vars <- rep(1., ncol(df))
        names(vars) <- predictor_cols
    }
    
    t0 <- printlog('Fitting and evaluating linear models.', time=T)
    if (is_binary){

        N <- nrow(df)
        N_CAS <- sum(df$Outcome == 1)
        N_CON <- sum(df$Outcome == 0)
        STUDY_PREV <- N_CAS / N

        # logistic link function
        models <- lapply(predictor_cols, function(x){
            df[,slim_glm(.SD),.SDcols=c(x,'Outcome')]
        })
        names(models) <- predictor_cols

        # confidence intervals 
        confint_coef_p95 <- t(sapply(models, function(x){confint.default(x)[2,]}))
        colnames(confint_coef_p95) <- c('BETA_CI_LOW', 'BETA_CI_HIGH')

        # get the summaries and delete large data
        model_summaries <- list()
        for (i in seq_along(models)){
            model_summaries[[i]] <- summary(models[[i]])
        }
        names(model_summaries) <- predictor_cols
        rm(models) # models is not needed anymore
        gc()

        # coefficients, standard errors and p-values
        coef_stderr_pval <- t(sapply(model_summaries, function(x){coef(x)[2,]}))
        colnames(coef_stderr_pval) <- c('BETA','SD','Z_VAL','PVAL')
        rm(model_summaries) # model_summaries is not needed anymore
        gc()
        t1 <- printlog('Completed fitting and processing glm models in ', time_taken(t0))

        # models on the observed scale
        models_obs <- lapply(predictor_cols, function(x){
            df[,slim_lm(.SD),.SDcols=c(x,'Outcome')]
        })
        names(models_obs) <- predictor_cols

        # Eq (7), Lee, S. H., M. E. Goddard, N. R. Wray and P. M. Visscher (2012)
        r2_obs <- as.matrix(sapply(models_obs, function(x){
            var(x$fitted.values) / (STUDY_PREV * (1 - STUDY_PREV))
        }))
        colnames(r2_obs) <- 'R2_OBS'

        # for binary phenotypes, confidence intervals on the observed scale (note: these might be inaccurate)
        confint_coef_p95_obs <- t(sapply(models_obs, function(x){confint.default(x)[2,]}))
        colnames(confint_coef_p95_obs) <- c('BETA_OBS_CI_LOW', 'BETA_OBS_CI_HIGH')

        # get the summaries and delete large data
        model_summaries_obs <- list()
        for (i in seq_along(models_obs)){
            model_summaries_obs[[i]] <- summary(models_obs[[i]])
            models_obs[[i]]$residuals <- NULL
            models_obs[[i]]$qr <- NULL
            models_obs[[i]]$effects <- NULL
        }
        names(model_summaries_obs) <- predictor_cols
        rm(models_obs) # models_obs not needed anymore
        gc()

        # for binary phenotypes, model coefficients, test statistics, and p-values on the observed scale
        coef_stderr_pval_obs <- t(sapply(model_summaries_obs, function(x){coef(x)[2,]}))
        colnames(coef_stderr_pval_obs) <- c('BETA_OBS','SD_OBS','T_VAL_OBS','PVAL_OBS')
        rm(model_summaries_obs) # models_summaries_obs not needed anymore
        gc()
        
        if (rescale_r2){
            r2_obs_ci <- confint_coef_p95_obs[,c('BETA_OBS_CI_LOW','BETA_OBS_CI_HIGH')]^2 / (STUDY_PREV * (1 - STUDY_PREV)) * cbind(vars, vars)
            colnames(r2_obs_ci) <- c('R2_OBS_CI_LOW', 'R2_OBS_CI_HIGH')
        } else {
            r2_obs_ci <- confint_coef_p95_obs[,c('BETA_OBS_CI_LOW','BETA_OBS_CI_HIGH')]^2 / (STUDY_PREV * (1 - STUDY_PREV))
            colnames(r2_obs_ci) <- c('R2_OBS_CI_LOW', 'R2_OBS_CI_HIGH')    
        }
        printlog('Completed fitting and processing lm models in ', time_taken(t1))

        if (opt$boot > 0){
            
            t0 <- printlog('Bootstrapping confidence intervals for R and R2 on the observed scale using ',opt$boot,' bootstrap samples.', time=T)
            
            # the bootstrap samples are not independent between scores
            boot_func <- function(dat, i){cor(dat[i,1],dat[i,2:ncol(dat)])}
            if (opt$threads > 1){
                boot_res <- boot(df[,c('Outcome', ..predictor_cols)], boot_func, R=opt$boot, strata=df$Outcome, ncpus = opt$threads, parallel = 'multicore')
            } else {
                boot_res <- boot(df[,c('Outcome', ..predictor_cols)], boot_func, R=opt$boot, strata=df$Outcome)    
            }
            
            r2_obs_ci_boot <- sapply(seq_along(predictor_cols), function(x){
                boot.ci(boot_res, type='perc', index = x, conf = 0.95)$percent
            })[4:5,]

            r2_obs_ci_boot<- t(r2_obs_ci_boot)
            row.names(r2_obs_ci_boot) <- predictor_cols
            
            r2_obs_ci_boot <- cbind(r2_obs_ci_boot, r2_obs_ci_boot^2)
            colnames(r2_obs_ci_boot) <- c('R_OBS_CI_BOOT_LOW', 'R_OBS_CI_BOOT_HIGH', 'R2_OBS_CI_BOOT_LOW', 'R2_OBS_CI_BOOT_HIGH')
            
            r2_obs_boot <- cbind(apply(boot_res$t, MARGIN = 2, median, na.rm=T),
                                 apply(boot_res$t, MARGIN = 2, mean, na.rm=T),
                                 apply(boot_res$t, MARGIN = 2, sd, na.rm=T),
                                 apply(boot_res$t, MARGIN = 2, function(x){median(x^2, na.rm=T)}),
                                 apply(boot_res$t, MARGIN = 2, function(x){mean(x^2, na.rm=T)}))
            r2_obs_boot <- cbind(t(boot_res$t0), r2_obs_boot)
            
            row.names(r2_obs_boot) <- predictor_cols
            colnames(r2_obs_boot) <- c('R_OBS_BOOT_T0','R_OBS_BOOT_MEDIAN','R_OBS_BOOT_MEAN','R_OBS_BOOT_SD','R2_OBS_BOOT_MEDIAN','R2_OBS_BOOT_MEAN')
            
            printlog('Boostrapping finished, runtime = ',time_taken(t0),time=T)
        }
        
        # AUC
        t0 <- printlog('Calculating AUC', time=T)
        auc_ci <- sapply(predictor_cols, function(x){
            suppressMessages(df[,ci.auc(formula = Outcome ~ . , data = data.frame(.SD), conf.level=0.95, method="delong"), .SDcols = c(x, 'Outcome')])
        })
        row.names(auc_ci) <- c('AUC_CI_LOW','AUC_MEDIAN','AUC_CI_HIGH')
        auc_ci <- t(auc_ci)
        printlog('Finished calculating AUC, runtime = ',time_taken(t0),time=T)

        
        if (opt$boot > 0){
            results <- Reduce(matrix_merge, list(coef_stderr_pval, confint_coef_p95, coef_stderr_pval_obs, confint_coef_p95_obs, r2_obs, r2_obs_ci, r2_obs_boot, r2_obs_ci_boot, auc_ci))
        } else {
            results <- Reduce(matrix_merge, list(coef_stderr_pval, confint_coef_p95, coef_stderr_pval_obs, confint_coef_p95_obs, r2_obs, r2_obs_ci, auc_ci))
        }
        
        # OR
        results[,OR:=exp(BETA)]
        results[,OR_CI_LOW:=exp(BETA_CI_LOW)]
        results[,OR_CI_HIGH:=exp(BETA_CI_HIGH)]
        
        # SAMPLE SIZE
        results$N <- N
        results$N_CAS <- N_CAS
        results$N_CON <- N_CON

    } else {

        # linear model (observed scale)
        models <- lapply(predictor_cols, function(x){
            df[,slim_lm(.SD),.SDcols=c(x,'Outcome')]
        })
        names(models) <- predictor_cols

        confint_coef_p95 <- t(sapply(models, function(x){confint.default(x)[2,]}))
        colnames(confint_coef_p95) <- c('BETA_CI_LOW', 'BETA_CI_HIGH')

        # get the summaries and delete large data
        model_summaries <- list()
        for (i in seq_along(models)){
            model_summaries[[i]] <- summary(models[[i]])
            models[[i]]$residuals <- NULL
            models[[i]]$qr <- NULL
            models[[i]]$effects <- NULL
        }
        names(model_summaries) <- predictor_cols
        rm(models) # models not needed anymore
        gc()

        # coefficients, standard errors and p-values
        coef_stderr_pval <- t(sapply(model_summaries, function(x){coef(x)[2,]}))
        colnames(coef_stderr_pval) <- c('BETA','SD','T_STAT','PVAL')

        # variance explained
        r2_obs <- as.matrix(sapply(model_summaries, function(x){
            x$r.squared
        }))
        colnames(r2_obs) <- c('R2_OBS')
        rm(model_summaries) # model_summaries not needed anymore
        gc()

        if (rescale_r2){
            r2_obs_ci <- confint_coef_p95[,c('BETA_CI_LOW','BETA_CI_HIGH')]^2 * cbind(vars, vars) / var(df$Outcome) 
            colnames(r2_obs_ci) <- c('R2_OBS_CI_LOW', 'R2_OBS_CI_HIGH')
        } else {
            r2_obs_ci <- confint_coef_p95[,c('BETA_CI_LOW','BETA_CI_HIGH')]^2 / var(df$Outcome)
            colnames(r2_obs_ci) <- c('R2_OBS_CI_LOW', 'R2_OBS_CI_HIGH')    
        }

        results <- Reduce(matrix_merge, list(coef_stderr_pval, confint_coef_p95, r2_obs, r2_obs_ci))    
        results$N <- nrow(df)
        
        printlog('Completed fitting and processing lm models in ', time_taken(t0))

    }
    # printlog('completed fitting and evaluating linear models, runtime = ',time_taken(t0), time=T)
    
    
    setnames(results,'rn','predictor')
    
    # by checking BETA, we can see if the confidence interval for the variance explained (R2) contains 0
    results[sign(BETA_CI_LOW) != sign(BETA_CI_HIGH), c('R2_OBS_CI_LOW','R2_OBS_CI_HIGH'):=list(0, pmax(R2_OBS_CI_HIGH, R2_OBS_CI_LOW))]
    
    if (opt$boot > 0 & is_binary){
        # by checking R, we can see if the confidence interval for the variance explained (R2) contains 0
        results[sign(R_OBS_CI_BOOT_LOW) != sign(R_OBS_CI_BOOT_HIGH), c('R2_OBS_CI_BOOT_LOW','R2_OBS_CI_BOOT_HIGH'):=list(0, pmax(R2_OBS_CI_BOOT_LOW, R2_OBS_CI_BOOT_HIGH))]
    }
    
    if (!is.na(opt$sex)){
        results$N_SEX0 <- sum(df$sex == 0)
        results$N_SEX1 <- sum(df$sex == 1)
    }
    
    return(results)
        
}

if (!is.na(opt$train_test)){
    t0 <- printlog('Calculating metrics within train-test splits.', time=T)
    results_train <- calculate_metrics(profiles[split=='train'])[,split:='train']
    t0 <- printlog('Finished calculating metrics for training set, runtime = ',time_taken(t0), time=T)
    results_test <- calculate_metrics(profiles[split=='test'], rescale_r2 = T)[,split:='test']
    printlog('Finished calculating metrics for test set, runtime = ',time_taken(t0),time=T)
    cat('\n')
    results_full <- rbindlist(list(results_train, results_test))
} else {
    t0 <- printlog('Calculating metrics for full dataset.', time=T)
    results_full <- calculate_metrics(profiles)[,split:='all']
    printlog('Finished calculating metrics for full dataset, runtime = ',time_taken(t0),'\n', time=T)
}

printlog('Writing metrics to ',paste0(opt$out,'.metrics.tsv'),'\n')
fwrite(results_full, paste0(opt$out,'.metrics.tsv'), sep='\t', na='NA', quote=F)

if (!is.na(opt$sex)){
    if (!all(profiles$sex == profiles$sex[1])){
        t0 <- printlog('Calculating metrics within sexes', time=T)
        results_sex0 <- calculate_metrics(profiles[sex==0], rescale_r2 = T)[,split:='sex0']
        t0 <- printlog('Finished calculating metrics for sex == 0, runtime = ',time_taken(t0), time=T)
        results_sex1 <- calculate_metrics(profiles[sex==1], rescale_r2 = T)[,split:='sex1']
        printlog('Finished calculating metrics for sex == 1, runtime = ',time_taken(t0), time=T)
        results_sex <- rbindlist(list(results_sex0, results_sex1))
        printlog('Writing sex-stratified metrics to ',paste0(opt$out,'.metrics_sex.tsv\n'))
        fwrite(results_sex, paste0(opt$out,'.metrics_sex.tsv'), sep='\t', na='NA', quote=F)
        cat('\n')
    } else {
        printlog('Warning: Sex is invariant. Skipping sex-specific analysis.\n')
    }
}

if (opt$calc_cor){
    t0 <- printlog('Calculating score-score correlations...',time=T)
} else {
    printlog('Analysis finished, runtime = ', time_taken(start.time), time=T)
    quit(save = 'no')
}

if (opt$sample_cor != -1){
    if (opt$sample_cor / nrow(profiles) > 1.){
        cat(paste0('--sample_cor has no effect because number of requested samples (',opt$sample_cor,') is larger than number of samples (',nrow(profiles),')\n'))
    } else {
        cat('Sampling to',opt$sample_cor,'individuals to calculate score-score correlations.\n')
        train_ind <- createDataPartition(profiles$Outcome, times = 1, list = F, p = opt$sample_cor/nrow(profiles))
        profiles <- profiles[train_ind]
    }
}

cors <- cor(profiles[,..predictor_cols])
printlog('Finished calculating score-score correlations, runtime = ',time_taken(t0))

printlog('Writing score-score correlations to ', paste0(opt$out,'.scor.tsv'),'\n')
write.table(cors, paste0(opt$out,'.scor.tsv'), sep='\t', na = 'NA', quote=F, row.names=T, col.names=NA)

printlog('Analysis finished, runtime = ', time_taken(start.time), time=T)



