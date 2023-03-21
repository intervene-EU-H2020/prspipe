suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))

infile <- commandArgs(trailingOnly = T)[1]

studies <- fread(infile)

# script to export a merged file containing the "automatic" methods
# takes one positional argument: the studies.tsv file

for (i in 1:nrow(studies)){
    
    study_id <- studies$study_id[i]
    
    all_scores <- data.frame(file=list.files('prs',pattern = '.*score.gz', full.names = T, recursive = T))
    all_scores$study <- sapply(strsplit(all_scores$file, '/'), '[[', 3)
    all_scores$method <- sapply(strsplit(all_scores$file, '/'), '[[', 2)
    all_scores <- filter(all_scores, study == study_id)

    heads <- sapply(all_scores$file, function(x){ fread(cmd = paste0(paste0('zcat ',x,' |  head -n 2')))})
    names(heads) <- all_scores$method

    selected_scores <- list()

    # lassosum

    lassosum_log_path <- list.files(paste0('prs/lassosum/',study_id,'/'), pattern = '\\.log$', full.names = T, recursive = F)

    if (length(lassosum_log_path) == 1){

        cat('Guessed lassosum log path : ',lassosum_log_path,'\n')
        cat('Reading best pseudovalidation parameters...\n')

        pseudoval_params<-system2(command='grep',args=c('-A','2','"Pseudovalidated parameters:"',lassosum_log_path),stdout = T, wait=T)

        stopifnot(startsWith(pseudoval_params[2],'s ='))
        stopifnot(startsWith(pseudoval_params[3],'lambda ='))

        lassosum_s <- gsub('s = ',replacement = 's',pseudoval_params[2],fixed = T)
        lassosum_lambda <- gsub('lambda = ',replacement = 'lambda',pseudoval_params[3],fixed = T)

    } else {
        cat('Could not guess lassosum log file path from input filename.\n')
        cat('Can\'t determine lassosum pseudovalidation model.\n')
    }

    selected_scores[['lassosum']] <- grep(paste0('SCORE_',lassosum_s,'_',gsub('[0-9]$','[0-9]+$',lassosum_lambda)), colnames(heads$lassosum))
    selected_scores[['lassosum']]

    # prscs
    selected_scores[['prscs']] <- grep('auto',colnames(heads$prscs))
    selected_scores[['prscs']]

    # sbayesr
    # only has a single score
    selected_scores[['sbayesr']] <- 3

    # dbslmm
    # only has a single score
    selected_scores[['dbslmm']] <- 3

    selected_scores[['ldpred2']] <- grep('auto|inf', colnames(heads$ldpred2))

    ##
    # MegaPRS
    ##

    megaprs_pseudoval_path <- list.files(paste0('prs/megaprs/',study_id,'/'), pattern = '\\.pseudo.+\\.txt$', full.names = T, recursive = F)

    if (length(megaprs_pseudoval_path)){

        cat('Guessed megaprs pseudovalidation log path : ',megaprs_pseudoval_path,'\n')
        cat('Reading best pseudovalidation parameters...\n')

        pseudoval_params<-fread(megaprs_pseudoval_path, header=F)

        pat <- paste0('ldak.Model', gsub('Score_','',pseudoval_params[which.max(V2)]$V1),'$')

    } else {
        cat('Warning: Can\'t determine MegaPRS pseudovalidation model.\n')
    }

    selected_scores[['megaprs']] <- grep(pat, colnames(heads$megaprs))


    all_scores <- filter(all_scores, method != 'pt_clump')

    scores <- lapply(1:nrow(all_scores), function(x){

        r <- all_scores[x,]

        i_select <- c(1,2,selected_scores[[r$method]])


        data <- fread(r$file, select = i_select)

        setnames(data, gsub('SCORE', r$method, colnames(data)))

        return(data)
    })

    ref <- fread('resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping_cached.tsv.gz', sep='\t', select=c('chr','pos_hg38','pos_hg19','rsid','a1','a2'))

    match_ref <- function(x){

        score_cols <- colnames(x)[3:ncol(x)]

        if (any(is.na(score_cols))){
            print(head(x))
            stop()
        }

        mrg <- merge(ref, x, by.x = 'rsid', by.y = 'SNP')

        for (s in score_cols){
            mrg[a2==A1,(s):=-1*get(s)]
        }

        ret_cols <- c('rsid','a1',score_cols)

        return(mrg[,..ret_cols])
    }

    scores <- lapply(scores, match_ref)

    scores_mrg <- Reduce(function(x,y){merge(x,y,by=c('rsid','a1'),all.x = T,all.y = T)}, scores)
    score_cols <- colnames(scores_mrg)[3:ncol(scores_mrg)]

    for (s in score_cols){
        scores_mrg[is.na(get(s)),(s):=0.]
    }

    scores_mrg <- merge(ref, scores_mrg, by=c('rsid','a1'))[order(chr, pos_hg38, a1)]

    setcolorder(scores_mrg, c('chr','pos_hg38','pos_hg19','rsid','a1','a2',score_cols))

    outfile <- paste0(study_id,'-auto_merged.tsv.gz')
    fwrite(scores_mrg, file = outfile, sep='\t')

    print(outfile)
}