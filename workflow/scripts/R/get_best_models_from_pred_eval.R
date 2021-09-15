#!/usr/bin/env RScript


suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(repr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))


option_list = list(
make_option("--pred_eval", action="store", default=NA, type='character',
        help="path to pred_eval.txt file"),
make_option("--drop", action="store", default=T, type='logical', 
        help="If FALSE, will report MultiPRS even if identical to Inf or Pseudoval models for SBLUP, DBSLMM")
)

opt = parse_args(OptionParser(option_list=option_list))
# predictors follow these patterns (could be expanded...)

predictors <- 'pt.clump|lassosum|prscs|sblup|sbayesr|dbslmm|ldpred2*'

# whether to drop MultiPRS from plots where it is (almost) identical to Inf/Auto
drop = TRUE

print(opt)

res_pred_eval <- fread(opt$pred_eval)

res_pred_eval$group<-gsub('\\.PredFile.*','',gsub('_group$','',res_pred_eval$Model))

res_pred_eval$method <- str_extract(res_pred_eval$Model, pattern = predictors)
res_pred_eval[is.na(method), method:='All']
res_pred_eval$is_multi <- 'no'
res_pred_eval[!grepl('File', Model), is_multi := 'yes']
res_pred_eval$select <- FALSE

is_categorical <- 'Cross_AUC' %in% colnames(res_pred_eval)

res_pred_eval$method_type <- ''
res_pred_eval$label <- ''

source('workflow/scripts/R/source_config.R')

##
# pT + clump
##
res_pred_eval[grepl(pattern = 'pt.clump.*PredFile', Model), c('label','method_type'):=list('pT+clump.10FCVal','CV')]
res_pred_eval[(method == 'pt.clump') & (is_multi == 'yes'), c('label','method_type'):=list('pT+clump.MultiPRS','MultiPRS')]
res_pred_eval[grepl(pattern = 'lassosum.*PredFile', Model), c('label','method_type'):=list('lassosum.10FCVal','CV')]
res_pred_eval[(method == 'lassosum') & (is_multi == 'yes'), c('label','method_type'):=list('lassosum.MultiPRS','MultiPRS')]

##
# lassosum
##
# skip lassosum pseudoval as it depends on the log file and can't be determined from the input files here

res_pred_eval[grepl(pattern = 'lassosum.*PredFile', Model), c('label','method_type'):=list('lassosum.10FCVal','CV')]
res_pred_eval[(method == 'lassosum') & (is_multi == 'yes'), c('label','method_type'):=list('lassosum.MultiPRS','MultiPRS')]

##
# PRScs
##
res_pred_eval[grepl(pattern = 'prscs.*PredFile.*phi[0-9]', Model), c('label','method_type'):=list('PRScs.10FCVal','CV')]
res_pred_eval[grepl(pattern = 'prscs.*PredFile.*phiauto', Model), c('label','method_type'):=list('PRScs.PseudoVal','PseudoVal')]
res_pred_eval[(method == 'prscs') & (is_multi == 'yes'), c('label','method_type'):=list('PRScs.MultiPRS','MultiPRS')]


##
# SBLUP
##
res_pred_eval[grepl(pattern = 'sblup.*PredFile.*SBLUP', Model), c('label','method_type'):=list('SBLUP.Inf','Inf')]

# since SBLUP only gives a single model, the "Multi-PRS" are identical to the single models
if (drop){
    res_pred_eval <-  res_pred_eval[!((method == 'sblup') & (is_multi == 'yes'))]
} else {
    res_pred_eval[(method == 'sblup') & (is_multi == 'yes'), c('label','method_type'):=list('SBLUP.MultiPRS*','MultiPRS') ] 
}

##
# SBayesR
##

res_pred_eval[grepl(pattern = 'sbayesr.*PredFile', Model), c('label','method_type'):=list('SBayesR.PseudoVal*','PseudoVal')] # SBayesR, even though it's a "pseudovalidated" method, was run with 2 sets of parameters
res_pred_eval[(method == 'sbayesr') & (is_multi == 'yes'), c('label','method_type'):=list('SBayesR.MultiPRS*','MultiPRS')] 

##
# LDPred1
##

res_pred_eval[grepl(pattern = 'ldpred.*PredFile.*p[0-9]', Model), c('label','method_type'):=list('LDPred.10FCVal','CV')]
res_pred_eval[grepl(pattern = 'ldpred.*PredFile.*inf', Model), c('label','method_type'):=list('LDPred.Inf','Inf')]
res_pred_eval[(method == 'ldpred') & (is_multi == 'yes'), c('label','method_type'):=list('LDPred.MultiPRS','MultiPRS')]

##
# LDPred2
##

res_pred_eval[grepl(pattern = 'ldpred2.*PredFile', Model) & !grepl(pattern = 'auto|inf', Model), c('label','method_type'):=list('LDPred2.10FCVal','CV')]
res_pred_eval[grepl(pattern = 'ldpred2.*PredFile.*inf', Model), c('label','method_type'):=list('LDPred2.Inf','Inf')]
# note: these can be missing due to infinite values!
res_pred_eval[grepl(pattern = 'ldpred2.*PredFile.*auto', Model), c('label','method_type'):=list('LDPred2.PseudoVal','PseudoVal')]

res_pred_eval[(method == 'ldpred2') & (is_multi == 'yes'), c('label','method_type'):=list('LDPred2.MultiPRS','MultiPRS')]

##
# DBSLMM
##

res_pred_eval[grepl(pattern = 'dbslmm.*PredFile', Model), c('label','method_type'):=list('DBSLMM.PseudoVal','PseudoVal')]

# since DBSLMM only gives a single model, the "Multi-PRS" are identical to the single models

if (drop){
    res_pred_eval <- res_pred_eval[!((method == 'dbslmm') & (is_multi == 'yes'))]
} else {
    res_pred_eval <- res_pred_eval[(method == 'dbslmm') & (is_multi == 'yes'), c('label','method_type'):=list('DBSLMM.MultiPRS*','MultiPRS')]
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

# export data to TSV
if (is_categorical){
    fwrite(res_eval[,list(group, Method, label, tag, method_type, is_multi, CrossVal_R, CrossVal_R_SE, CrossVal_pval, CrossVal_N, IndepVal_R, IndepVal_R_SE, IndepVal_pval, IndepVal_N, CrossVal_OR, CrossVal_LowCI, CrossVal_HighCI, Cross_LiabR2, Cross_AUC, CrossVal_Ncas, CrossVal_Ncon, IndepVal_OR, IndepVal_LowCI, IndepVal_HighCI, Indep_LiabR2, Indep_AUC, IndepVal_Ncas, IndepVal_Ncon)], paste0(outpath, '/best_models.tsv'), sep='\t')
} else {
    fwrite(res_eval[,list(group, Method, label, tag, method_type, is_multi, CrossVal_R, CrossVal_R_SE, CrossVal_pval, CrossVal_N, IndepVal_R, IndepVal_R_SE, IndepVal_pval, IndepVal_N)], paste0(outpath, '/best_models.tsv'), sep='\t')
}
