suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))

option_list = list(
make_option("--pred_eval", action="store", default=NA, type='character',
        help="path to pred_eval.txt file"),
make_option("--drop", action="store", default=T, type='logical', 
        help="If FALSE, will report MultiPRS even if identical to Inf or Pseudoval models for DBSLMM, SBayesR")
)

predictors <- 'pt.clump|lassosum|prscs|sblup|sbayesr|dbslmm|ldpred2|megaprs|All'

opt = parse_args(OptionParser(option_list=option_list))


# whether to drop MultiPRS from plots where it is (almost) identical to Inf/Auto
drop = opt$drop

print(opt)
res_pred_eval <- fread(opt$pred_eval)


res_pred_eval$group<-gsub('\\.PredFile.*','',gsub('_group$','',res_pred_eval$Model))

res_pred_eval$method <- str_extract(res_pred_eval$Model, pattern = predictors)

res_pred_eval$is_multi <- 'no'
res_pred_eval[!grepl('[Ff]ile', Model), is_multi := 'yes']
res_pred_eval$select <- FALSE

is_categorical <- 'Cross_AUC' %in% colnames(res_pred_eval)

res_pred_eval$method_type <- ''
res_pred_eval$label <- ''

source('workflow/scripts/R/source_config.R')

##
# pT + clump
##
res_pred_eval[grepl(pattern = 'pt.clump.*PredFile', Model), c('label','method_type'):=list('pT+clump.CV','CV')]
res_pred_eval[(method == 'pt.clump') & (is_multi == 'yes'), c('label','method_type'):=list('pT+clump.MultiPRS','MultiPRS')]


##
# lassosum
##

res_pred_eval[grepl(pattern = 'lassosum.*PredFile', Model), c('label','method_type'):=list('lassosum.CV','CV')]
res_pred_eval[(method == 'lassosum') & (is_multi == 'yes'), c('label','method_type'):=list('lassosum.MultiPRS','MultiPRS')]

pheno <- strsplit(str_extract(opt$pred_eval, pattern = 'PRS_evaluation/([_\\-\\.]|[:alnum:])+'),'/')[[1]][2]

if (!is.na(pheno)){
    
    lassosum_log_path <- list.files(paste0('prs/lassosum/',pheno,'/'), pattern = '\\.log$', full.names = T, recursive = F)
    
    if (length(lassosum_log_path) == 1){
    
        cat('Guessed lassosum log path : ',lassosum_log_path,'\n')
        cat('Reading best pseudovalidation parameters...\n')
        
        pseudoval_params<-system2(command='grep',args=c('-A','2','"Pseudovalidated parameters:"',lassosum_log_path),stdout = T, wait=T)
        
        stopifnot(startsWith(pseudoval_params[2],'s ='))
        stopifnot(startsWith(pseudoval_params[3],'lambda ='))
        
        lassosum_s <- gsub('s = ',replacement = 's',pseudoval_params[2],fixed = T)
        lassosum_lambda <- gsub('lambda = ',replacement = 'lambda',pseudoval_params[3],fixed = T)
        
        if (nchar(lassosum_lambda) > 11){
            lassosum_lambda <- gsub('[0-9]$','[0-9]+_group',lassosum_lambda) 
        } else {
            lassosum_lambda <- paste0(lassosum_lambda, '_group')
        }
        
        pat <- paste0(lassosum_s,'.',lassosum_lambda)
        row <- res_pred_eval[grepl(pattern = pat, Model) & grepl(pattern = 'lassosum.*PredFile', Model)]

        print(pat)
        print(res_pred_eval[grepl(pattern = 'lassosum.*PredFile', Model)])

        if (nrow(row) == 1){    
            row[,c('label','method_type'):=list('lassosum.PseudoVal','PseudoVal')]
            res_pred_eval <- rbindlist(list(res_pred_eval, row))
        } else {
            cat('Can\'t determine lassosum pseudovalidation model.\n')
            stop()
        }
        rm(pat, row, pseudoval_params)
    }
} else {
    cat('Could not guess lassosum log file path from input filename.\n')
    cat('Can\'t determine lassosum pseudovalidation model.\n')
}


##
# PRScs
##

res_pred_eval[grepl(pattern = 'prscs.*PredFile.*', Model), c('label','method_type'):=list('PRScs.CV','CV')]
res_pred_eval[(method == 'prscs') & (is_multi == 'yes'), c('label','method_type'):=list('PRScs.MultiPRS','MultiPRS')]

# pseudovalidation 
row <- res_pred_eval[grepl(pattern = 'prscs.*PredFile.*phi\\.auto', Model)]
if (nrow(row) == 1){
    row[,c('label','method_type'):=list('PRScs.PseudoVal','PseudoVal')]
    res_pred_eval <- rbindlist(list(res_pred_eval, row))
    rm(row)
} else {
    cat('Can\'t determine PRScs pseudovalidation model.\n')
}


##
# SBayesR
##

res_pred_eval[grepl(pattern = 'sbayesr.*PredFile', Model), c('label','method_type'):=list('SBayesR.PseudoVal','PseudoVal')]
res_pred_eval[(method == 'sbayesr') & (is_multi == 'yes'), c('label','method_type'):=list('SBayesR.MultiPRS','MultiPRS')] 


#
# LDPred2
##

res_pred_eval[grepl(pattern = 'ldpred2.*PredFile', Model), c('label','method_type'):=list('LDPred2.CV','CV')]
res_pred_eval[(method == 'ldpred2') & (is_multi == 'yes'), c('label','method_type'):=list('LDPred2.MultiPRS','MultiPRS')]


# note: these can be missing due to infinite values!
row <- res_pred_eval[grepl(pattern = 'ldpred2.*PredFile.*inf', Model)]
if (nrow(row) == 1){
    row[, c('label','method_type'):=list('LDPred2.Inf','Inf')]
    res_pred_eval <- rbindlist(list(res_pred_eval, row))
} else {
    cat('Can\'t determine LDPred2 Inf model.\n')
}
rm(row)

row <- res_pred_eval[grepl(pattern = 'ldpred2.*auto', Model),]
if (nrow(row) == 1){
    row[, c('label','method_type'):=list('LDPred2.PseudoVal','PseudoVal')]
    res_pred_eval <- rbindlist(list(res_pred_eval, row))
} else {
    cat('Can\'t determine LDPred2 pseudovalidation model.\n')
}
rm(row)


##
# DBSLMM
##

res_pred_eval[grepl(pattern = 'dbslmm.*PredFile', Model), c('label','method_type'):=list('DBSLMM.PseudoVal','PseudoVal')]

# since DBSLMM only gives a single model, the "Multi-PRS" are identical to the single models

if (drop){
    res_pred_eval <- res_pred_eval[!((method == 'dbslmm') & (is_multi == 'yes'))]
} else {
    res_pred_eval <- res_pred_eval[(method == 'dbslmm') & (is_multi == 'yes'), c('label','method_type'):=list('DBSLMM.MultiPRS','MultiPRS')]
}


##
# MegaPRS
##

res_pred_eval[grepl(pattern = 'megaprs.*PredFile', Model), c('label', 'method_type'):=list('MegaPRS.CV', 'CV')]
res_pred_eval[(method == 'megaprs') & (is_multi == 'yes'), c('label', 'method_type'):=list('MegaPRS.MultiPRS', 'MultiPRS')]

if (!is.na(pheno)){
    
    megaprs_pseudoval_path <- list.files(paste0('prs/megaprs/',pheno,'/'), pattern = '\\.pseudo.+\\.txt$', full.names = T, recursive = F)
    
    if (length(megaprs_pseudoval_path)){
    
        cat('Guessed megaprs pseudovalidation log path : ',megaprs_pseudoval_path,'\n')
        cat('Reading best pseudovalidation parameters...\n')
        
        pseudoval_params<-fread(megaprs_pseudoval_path, header=F)
        
        pat <- paste0('ldak.Model', gsub('Score_','',pseudoval_params[which.max(V2)]$V1))
        row <- res_pred_eval[grepl(pattern = pat, Model) & grepl(pattern = 'megaprs*', Model)]
        
        if (nrow(row) == 1){    
            row[,c('label','method_type'):=list('MegaPRS.PseudoVal','PseudoVal')]
            res_pred_eval <- rbindlist(list(res_pred_eval, row))
        } else {
            cat('Can\'t determine megaprs pseudovalidation model.\n')
        }
        rm(pat, row, pseudoval_params)
    }
} else {
    cat('Could not guess lassosum log file path from input filename.\n')
    cat('Can\'t determine lassosum pseudovalidation model.\n')
}

res_pred_eval[method=='All' ,  c('label','method_type'):=list('All.MultiPRS','MultiPRS')]
res_pred_eval[,select:=CrossVal_R == max(CrossVal_R),by=list(group,method_type)]

res_eval <- res_pred_eval[(select),]
res_eval$Method <- res_eval$method
res_eval$tag <- res_eval$Model
res_eval$Model<-factor(res_eval$method_type, level=c('MultiPRS','CV','PseudoVal','Inf','All'))
res_eval$Test <- res_eval$label

outpath <- paste0(dirname(opt$pred_eval))


cat('writing output to', paste0(outpath, '/best_models.tsv'), '\n')

res_eval[,select:=NULL]
res_eval[,Test:=NULL]

fwrite(res_eval,  paste0(outpath, '/best_models.tsv'), sep='\t', row.names = F, quote=F)
