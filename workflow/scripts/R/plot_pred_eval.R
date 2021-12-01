#!/usr/bin/env RScript


suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))


option_list = list(
make_option("--pred_eval", action="store", default=NA, type='character',
        help="path to pred_eval.txt file"),
make_option("--pheno_name", action="store", default=NA, type='character',
        help="name of phenotype to label plots"),
make_option("--drop", action="store", default=T, type='logical', 
        help="If FALSE, will report MultiPRS even if identical to Inf or Pseudoval models for SBLUP, DBSLMM")
)

opt = parse_args(OptionParser(option_list=option_list), args=c('--pred_eval', './results/estonia/GCST002783/GCST002783.AllMethodComp.pred_eval.txt'))
# predictors follow these patterns (could be expanded...)

out_prefix <- paste0(gsub('\\.AllMethodComp\\.pred_eval\\.txt$', '',opt$pred_eval))
study<-basename(out_prefix)

if (is.na(opt$pheno_name)){
    opt$pheno_name <- study
    cat(paste0('No phenotype name defined, defaulting to study name (',study,')\n'))
}

cat(paste0('Plotting performances for study ',study),'\n')


predictors <- 'pt.clump|lassosum|prscs|sblup|sbayesr|dbslmm|ldpred2*'

# whether to drop MultiPRS from plots where it is (almost) identical to Inf/Auto
drop = TRUE

cat('options in place:\n')
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

# res_eval # these are the best models...

##
# parameter ranges as defined in the separate polygenic score file creator scripts
##

# pT + clump
# default values for the p-value thresholds
# '1e-8,1e-6,1e-4,1e-2,0.1,0.2,0.3,0.4,0.5,1'

# p-value threshold is the only variable that is optimized 

ptclump_i <- grep('pt\\.clump\\.PredFile[0-9]+', res_pred_eval$Model, value = F)
ptclump_predfiles <- res_pred_eval$Model[ptclump_i]

ptclump_pTs <- gsub('pt\\.clump\\.PredFile[0-9]+\\.[A-Z,a-z,0-9]+\\.', '', ptclump_predfiles)
ptclump_pTs <- gsub('_group','',ptclump_pTs)
ptclump_pTs <- gsub('e\\.','e-',ptclump_pTs)
ptclump_pTs <- as.numeric(ptclump_pTs)

options(repr.plot.width=4, repr.plot.height=3, repr.plot.res=300)

plotdata <- res_pred_eval[ptclump_i]
plotdata$pvalue_threshold <- ptclump_pTs

# plotdata<-melt(plotdata, id.vars = c('pheno','study','pvalue_threshold'), measure.vars = c('CrossVal_R', 'IndepVal_R'))
# ggplot(plotdata, aes(x=log10(pvalue_threshold), y=value, col=variable))  + geom_line(aes(group=variable)) + geom_point() + facet_wrap(~pheno, nrow = 1) + labs(title = 'pT + clump') + xlab('log10(p-value threshold)') + ylab('R')
plotdata<-melt(plotdata, id.vars = c('Model','pvalue_threshold'), measure.vars = c('CrossVal_R', 'IndepVal_R'))

pdf(paste0(out_prefix, '_performance_plots.pdf'),width = 4.5, 3)
ggplot(plotdata, aes(x=log10(pvalue_threshold), y=value, col=variable))  + geom_line(aes(group=variable)) + geom_point() + labs(title = 'pT + clump') + xlab('log10(p-value threshold)') + ylab('R')

plotdata$method <- 'pT+clump'
plotdata$study <- study

fwrite(plotdata, file = paste0(out_prefix, '.pdata_ptclump.txt'))

# lassosum 
# infer from file names

# s -> mixing parameter (?)
# lambda -> weight of L1 penalty (?)

options(repr.plot.width=4, repr.plot.height=3, repr.plot.res=300)


lassosum_i <- grep('lassosum.PredFile', res_pred_eval$Model)
lassosum_lambda <- res_pred_eval$Model[lassosum_i]

lassosum_s <- as.numeric(gsub('\\.$','',gsub('s','',str_extract_all('s[0-2]\\.[0-9]*', string = lassosum_lambda))))
lassosum_lambda <- as.numeric(gsub('e\\.','e-',gsub('_group','',gsub('lambda','',str_extract_all('lambda.*$', string = lassosum_lambda)))))

plotdata <- res_pred_eval[lassosum_i]
plotdata$s <- lassosum_s
plotdata$lambda <- lassosum_lambda
# plotdata <- melt(plotdata,  measure.vars = c('CrossVal_R', 'IndepVal_R'), id.vars = c('pheno','study','s','lambda'))
plotdata <- melt(plotdata,  measure.vars = c('CrossVal_R', 'IndepVal_R'), id.vars = c('Model','s','lambda'))

# ggplot(plotdata, aes(y=s, x=log10(lambda)))+geom_point(aes(col=value, size=value)) + facet_wrap(~pheno, nrow = 2)

plotdata$s <- factor(plotdata$s, levels=sort(unique(plotdata$s)))
# ggplot(filter(plotdata, variable=='CrossVal_R'), aes(y=value, x=log10(lambda)))+geom_line(aes(group=s, col=s))+geom_point(aes(col=s)) + facet_wrap(~pheno, nrow = 2) + labs(title = 'lassosum - CV R') + xlab('log10(lambda)') + ylab('R')
ggplot(filter(plotdata, variable=='CrossVal_R'), aes(y=value, x=log10(lambda)))+geom_line(aes(group=s, col=s))+geom_point(aes(col=s)) + labs(title = 'lassosum - CV R') + xlab('log10(lambda)') + ylab('R')

ggplot(filter(plotdata, variable=='IndepVal_R'), aes(y=value, x=log10(lambda)))+geom_line(aes(group=s, col=s))+geom_point(aes(col=s)) + labs(title = 'lassosum - Indep R') + xlab('log10(lambda)') + ylab('R')

plotdata$method <- 'lassosum'
plotdata$study <- study

fwrite(plotdata, file = paste0(out_prefix, '.pdata_lassosum.txt'))

# PRScs 
# infer from file names
prscs_i <- grep('prscs\\.PredFile.*phi[0-9]+', res_pred_eval$Model)
prscs_phi <- res_pred_eval$Model[prscs_i]

prscs_phi <- gsub('.*phi','',prscs_phi)
prscs_phi <- gsub('_group','',prscs_phi)
prscs_phi <- gsub('e\\.', 'e-',prscs_phi)
prscs_phi <- as.numeric(prscs_phi)

plotdata <- res_pred_eval[prscs_i]
plotdata$phi <- prscs_phi

# plotdata <- melt(plotdata,  measure.vars = c('CrossVal_R', 'IndepVal_R'), id.vars = c('pheno','study','phi'))
plotdata <- melt(plotdata,  measure.vars = c('CrossVal_R', 'IndepVal_R'), id.vars = c('Model','phi'))

options(repr.plot.width=4, repr.plot.height=3, repr.plot.res=300)
# ggplot(plotdata, aes(x=log10(phi), y=value, col=variable))  + geom_line(aes(group=variable)) + geom_point() + facet_wrap(~pheno, nrow = 1) + labs(title = 'PRScs') + xlab('log10(phi)') + ylab('R')
ggplot(plotdata, aes(x=log10(phi), y=value, col=variable))  + geom_line(aes(group=variable)) + geom_point() + labs(title = 'PRScs') + xlab('log10(phi)') + ylab('R')

plotdata$method <- 'PRScs'
plotdata$study <- study

fwrite(plotdata, file = paste0(out_prefix, '.pdata_prscs.txt'))

# LDpred 
# infer from file names
options(repr.plot.width=4, repr.plot.height=3, repr.plot.res=300)


ldpred_i <- grep('ldpred\\.PredFile.*p[0-9]+', res_pred_eval$Model)
ldpred_p <- res_pred_eval$Model[ldpred_i]


ldpred_p <- gsub('.*LDpred\\.p','',ldpred_p)
ldpred_p <- gsub('_group','',ldpred_p)
ldpred_p <- gsub('e\\.', 'e-',ldpred_p)
ldpred_p <- as.numeric(ldpred_p)

plotdata <- res_pred_eval[ldpred_i]
plotdata$p <- ldpred_p
# plotdata <- melt(plotdata,  measure.vars = c('CrossVal_R', 'IndepVal_R'), id.vars = c('pheno','study','p'))
 plotdata <- melt(plotdata,  measure.vars = c('CrossVal_R', 'IndepVal_R'), id.vars = c('Model','p'))

# ggplot(plotdata, aes(x=log10(p), y=value, col=variable))  + geom_line(aes(group=variable)) + geom_point() + facet_wrap(~pheno, nrow = 1) + labs(title = 'LDpred') + xlab('log10(p)') + ylab('R')
ggplot(plotdata, aes(x=log10(p), y=value, col=variable))  + geom_line(aes(group=variable)) + geom_point() + labs(title = 'LDpred') + xlab('log10(p)') + ylab('R')

plotdata$method <- 'LDPred1'
plotdata$study <- study

fwrite(plotdata, file = paste0(out_prefix, '.pdata_ldpred.txt'))

# LDpred2 
# infer from file names
ldpred2_i <- which(grepl('ldpred2\\.PredFile.*', res_pred_eval$Model) & !(grepl('ldpred2\\.PredFile.*beta', res_pred_eval$Model)))
ldpred2_p <- res_pred_eval$Model[ldpred2_i]

# for (s in unique(md$study_id)){
for (s in unique(c('GCST002783'))){
    ldpred2_p <- gsub(s,'',ldpred2_p)
}

ldpred2_p <- gsub('.*PredFile[0-9]+\\.\\.','',ldpred2_p)
ldpred2_p <- gsub('_group','',ldpred2_p)

ldpred2_all <- ldpred2_p # keep full string stored
ldpred2_all <- gsub('e\\.','e-',ldpred2_all)

ldpred2_p <- str_extract(ldpred2_all, pattern = '^[0-9]+\\.*[0-9]*(e\\-)?[0-9]*')

ldpred2_h2 <- gsub(pattern = '(^[0-9](\\.[0-9]*)?e-[0-9]+\\.)|(^[0-9]+\\.[0-9]+\\.)','',ldpred2_all)
ldpred2_h2 <- gsub(pattern = '\\.sparse|\\.nosparse', '', ldpred2_h2)
ldpred2_h2 <- paste0(ifelse(grepl('\\.',ldpred2_h2), yes = '',no='0.'), ldpred2_h2)

ldpred2_sparse <- str_extract('sparse|nosparse', string = ldpred2_all)

plotdata <- res_pred_eval[ldpred2_i]
plotdata$h2 <- as.numeric(ldpred2_h2)
plotdata$p <- as.numeric(ldpred2_p)
plotdata$sparse <- ldpred2_sparse

# plotdata <- melt(plotdata,  measure.vars = c('CrossVal_R', 'IndepVal_R'), id.vars = c('pheno','study','p','h2','sparse'))
plotdata <- melt(plotdata,  measure.vars = c('CrossVal_R', 'IndepVal_R'), id.vars = c('Model','p','h2','sparse'))

options(repr.plot.width=4, repr.plot.height=3, repr.plot.res=300)
# ggplot(filter(plotdata, variable=='CrossVal_R', sparse == 'sparse'), aes(y=value, x=log10(p)))+geom_line(aes(group=h2, col=h2))+geom_point(aes(col=h2)) + facet_wrap(~pheno, nrow = 2) + labs(title = 'LDpred2 - CV R (sparse)') + xlab('log10(p)') + ylab('R')

ggplot(filter(plotdata, variable=='CrossVal_R', sparse == 'nosparse'), aes(y=value, x=log10(p)))+geom_line(aes(group=h2, col=h2))+geom_point(aes(col=h2)) + labs(title = 'LDpred2 - CV R (dense)') + xlab('log10(p)') + ylab('R')
ggplot(filter(plotdata, variable=='IndepVal_R', sparse == 'nosparse'), aes(y=value, x=log10(p)))+geom_line(aes(group=h2, col=h2))+geom_point(aes(col=h2)) + labs(title = 'LDpred2 - Indep R (dense)') + xlab('log10(p)') + ylab('R')


plotdata$method <- 'LDPred2'
plotdata$study <- study

fwrite(plotdata, file = paste0(out_prefix, '.pdata_ldpred2.txt'))

res_eval <- res_pred_eval[(select),]

method_mapping <- list(
    'pt.clump' = 'pT+clump',
    'lassosum' = 'lassosum',
    'All' = 'All',
    'prscs' = 'PRScs',
    'sblup' = 'SBLUP',
    'sbayesr' = 'SBayesR',
    'ldpred' = 'LDPred1',
    'ldpred2' = 'LDPred2',
    'dbslmm' = 'DBSLMM'
)

res_eval$Model_1 <- res_eval$Model # backup

res_eval$Method <- sapply(res_eval$method, function(x){method_mapping[[x]]})
res_eval$Method <- factor(res_eval$Method, levels=unique(res_eval$Method))
res_eval$Model<-factor(res_eval$method_type, level=c('MultiPRS','CV','PseudoVal','Inf'))
res_eval$Test <- res_eval$label

library(cowplot)

res_eval$pheno <- opt$pheno_name
md<-data.frame(name=c(opt$pheno_name)) # this is done weirdly here because I just copied code from a different script...

phenos <- unique(md$name)

options(repr.plot.width=6, repr.plot.height=4, repr.plot.res=300)


for (i in seq_along(phenos)){
    
    plotdata <- res_eval[pheno == phenos[i]]
    
    study <- plotdata$study[1]
    
    plotdata <- plotdata[,list(Method, Model, CrossVal_R, IndepVal_R, CrossVal_R_SE, IndepVal_R_SE)]
    
    for(j in unique(plotdata$Method)){
      for(k in unique(plotdata$Model)){
        if(dim(plotdata[Method == j & Model == k,])[1] > 0){
          next
        } else {
          dud<-data.frame('Method'=j,
                          'Model'=k,
                          'CrossVal_R'=NA,
                          'CrossVal_R_SE'=NA,
                          'IndepVal_R'=NA,
                          'IndepVal_R_SE'=NA
                          )
          plotdata<-rbind(plotdata,dud)
        }
      }
    }
    
    # export data to TSV
    # fwrite(plotdata, paste0('./results/UKBB/PRS_for_comparison/evaluation/',study,'/Association_withPRS/best_models.tsv'), sep='\t')
    
    
    # crossvalidation R
    print(ggplot(plotdata, aes(x=Method, y=CrossVal_R, fill=Model)) +
                          geom_bar(stat="identity", position = position_dodge(preserve = 'single'), width = 0.7) +
                          geom_errorbar(aes(ymin=CrossVal_R-CrossVal_R_SE, ymax=CrossVal_R+CrossVal_R_SE), width=.2, position=position_dodge(preserve = "single", width = 0.7)) +
                          geom_hline(yintercept = plotdata[(Method == 'pT+clump') & (Model == 'MultiPRS')]$CrossVal_R, linetype=2, alpha=0.5) +
                          coord_cartesian(ylim=c(min(plotdata$CrossVal_R-0.02), max(plotdata$CrossVal_R+0.025)), clip="on") +
                          labs(y="Correlation (SE)", x='', title=phenos[i]) +
                          theme_half_open() +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                          background_grid(major = 'y', minor = 'y') +
                          geom_vline(xintercept = 1.5, linetype="dotted") +
                          geom_vline(xintercept = 2.5, linetype="dotted") + 
                          geom_vline(xintercept = 3.5, linetype="dotted") + 
                          geom_vline(xintercept = 4.5, linetype="dotted") + 
                          geom_vline(xintercept = 5.5, linetype="dotted") +
                          geom_vline(xintercept = 6.5, linetype="dotted") +
                          geom_vline(xintercept = 7.5, linetype="dotted") +
                          geom_vline(xintercept = 8.5, linetype="dotted") +
                          theme(legend.position="top", legend.title = element_blank(), legend.box="vertical", legend.margin=margin()) +
                          guides(fill=guide_legend(nrow=2))) 

    
    # independent validation set R
    print(ggplot(plotdata,aes(x=Method, y=IndepVal_R, fill=Model)) +
                          geom_bar(stat="identity", position=position_dodge(preserve = "single"), width = 0.7) +
                          geom_errorbar(aes(ymin=IndepVal_R-IndepVal_R_SE, ymax=IndepVal_R+IndepVal_R_SE), width=.2, position=position_dodge(width = 0.7, preserve = "single")) +
                          geom_hline(yintercept = plotdata[(Method == 'pT+clump') & (Model == 'MultiPRS')]$IndepVal_R, linetype=2, alpha=0.5) +
                          labs(y="Correlation (SE)", x='', title=phenos[i]) +
                          coord_cartesian(ylim=c(min(plotdata$IndepVal_R-0.02), max(plotdata$IndepVal_R+0.025)), clip="on") +                          
                          theme_half_open() +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                          background_grid(major = 'y', minor = 'y') +
                          background_grid(major = 'y', minor = 'y') +
                          geom_vline(xintercept = 1.5, linetype="dotted") +
                          geom_vline(xintercept = 2.5, linetype="dotted") + 
                          geom_vline(xintercept = 3.5, linetype="dotted") + 
                          geom_vline(xintercept = 4.5, linetype="dotted") + 
                          geom_vline(xintercept = 5.5, linetype="dotted") +
                          geom_vline(xintercept = 6.5, linetype="dotted") +
                          geom_vline(xintercept = 7.5, linetype="dotted") +
                          geom_vline(xintercept = 8.5, linetype="dotted") +
                          theme(legend.position="top", legend.title = element_blank(), legend.box="vertical", legend.margin=margin()) +
                          guides(fill=guide_legend(nrow=2)))

}

dev.off()