#!/usr/bin/env Rscript


start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--model_rds", action="store", default=NA, type='character',
            help="Path to saved GenoPred glmnet-model"),
make_option("--hm3_1kg_mapping", action="store", default='resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping.tsv.gz', type='character'),
make_option("--prs_dir", action="store", default="prs/", type="character")
)

opt = parse_args(OptionParser(option_list=option_list))

if (!endsWith(opt$prs_dir, '/')){
    opt$prs_dir <- paste0(opt$prs_dir,'/')    
}

stopifnot(file.exists(opt$hm3_1kg_mapping))
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

predictor_groups_file <- paste0(dirname(model_rds),'/',study_id,'.AllMethodComp.predictor_groups')
cat(paste0('Using predictor groups file: ', predictor_groups_file, '\n'))
predictor_groups <- fread(predictor_groups_file)

score_files <- lapply(predictor_groups$group, function(x){list.files(paste0(opt$prs_dir,x),pattern = paste0('1KGPhase3.w_hm3.',study_id,'.score.gz'), full.names = T, recursive = T)})

for (i in seq_along(score_files)){
    if (length(score_files[[i]]) == 0){
        stop(paste0('Could not identify score file for group ', predictor_groups$group[i]))
    } else if (length(score_files[[i]]) > 1){
        stop(paste0('Found more than one score file for group ', predictor_groups$group[i], ' -> ', paste(score_files[[i]], collapse=';')))
    }
    
    cat(paste0('Found score file ', paste(score_files[[i]], collapse=','),' for group "',predictor_groups$group[i],'"\n'))

}

scale_files <- sapply(score_files, function(x){gsub('\\.score(\\.gz){0,1}$', paste0('.',ancestry,'.scale'), x)})

if (!all(file.exists(scale_files))){
    stop('Missing scale files.')    
}


predictor_groups$group <- gsub("[[:punct:]]", ".",predictor_groups$group)

score_files <- unlist(score_files)
scales <- lapply(scale_files, fread)

scores_by_group <- lapply(seq_along(score_files), function(i){
    
    group <- predictor_groups$group[i]
    
    selected <- grep(paste0('^Group_',group), names(non_zero_coef), value=TRUE)
    
    if (length(selected) == 0){
        cat('No scores selected for group',group)
        cat(', skipping...\n')
        return(NULL)
    } else {
        cat('Selected', length(selected), 'scores for group',group,'\n')    
    }
    
    predfilenum <- str_extract(string=selected, 'PredFile[0-9]+')
    
    if (!all(predfilenum == predfilenum[1])){
        # unable to determine the "PredFile" number
        print(predfilenum)
        stop(paste0('Could not determine PredFile for group ', group))    
    }
    
    predfilenum <- predfilenum[1]
    
    score_head <- fread(score_files[i], nrows = 1, header = T)
    score_columns <- colnames(score_head)
    
    new_names <- gsub('^SCORE',study_id,score_columns)
    new_names <- gsub("[-+]", ".", new_names)
    new_names <- paste0('Group_', group, '.', predfilenum, '.', new_names)

    if (!all(selected %in% new_names)){
        stop(paste0('Missing scores for group ', group))
    }
            
    matched_idx <- match(selected, new_names)
    read_cols <- c(1,2,matched_idx) # read only SNP, A1, and selected columns
    
    original_col_names <- score_columns[matched_idx]
    
    if (!all(original_col_names %in% scales[[i]]$Param)){
        # .scale file does not contain the score
        # print(scales[[i]]$Param)
        # print(original_col_names)
        stop(paste0('Some scores are missing from the .scale file for group "', group,'"'))    
    }

    score <- fread(score_files[i], header = T, select = read_cols)

    # rescale the effect-sizes by the 1KG standard deviation (as was done for the polygenic scores before fitting the elastic net model)
    SDs <- scales[[i]]$SD[match(original_col_names, scales[[i]]$Param)]
    
    # scaled_polygenic_scorer substitutes 0.000 SD values with 1e-6, we also reverse this substitution here
    if (any(SDs < 1e-6)){
        cat('Warning: some of the selected scores have very small reference standard deviations.\n')
        SDs <- ifelse(SDs < 1e-6, 1e-6, SDs)
    }
    
    for (c in seq_along(SDs)){
        col <- colnames(score)[c+2]
        score[,(col):=get(col)/SDs[c]]
    }
    
    snp_a1 <- score[,list(SNP,A1)]
    
    score[,c('SNP','A1'):=NULL]
    
    # multiply the scaled scores with their parameteres in the elastic net model, and sum up
    snp_a1$scoresum <- as.matrix(score) %*% non_zero_coef[selected]
    setnames(snp_a1, 'scoresum', paste0('Group_',group))
    
    return(snp_a1)
})

# recode the scores to use the same A1 alleles
hm3_1kg_mapping <- fread(opt$hm3_1kg_mapping, select = c('rsid','a1','a2'), header=T)

setnames(hm3_1kg_mapping,'rsid','SNP')

for (i in seq_along(scores_by_group)){
    
    if (length(scores_by_group[[i]]) == 0){
        next
    }
    
    group <- colnames(scores_by_group[[i]])[3]
    
    matching <- match(scores_by_group[[i]]$SNP, hm3_1kg_mapping$SNP)
    
    if (any(is.na(matching))){
        n_missing <- sum(is.na(matching))
        cat(paste0('Warning, ',n_missing,' SNPs are missing in the HapMap3-1KG reference. (',100*n_missing/nrow(scores_by_group[[i]]),'%)\n'))   
    }
    
    scores_by_group[[i]]$a1 <- hm3_1kg_mapping$a1[matching]
    scores_by_group[[i]]$a2 <- hm3_1kg_mapping$a2[matching]
    
    scores_by_group[[i]][a1 == A1, effect:=get(group)]
    scores_by_group[[i]][a2 == A1, effect:=-1*get(group)]
    
    scores_by_group[[i]][[group]] <- NULL
    setnames(scores_by_group[[i]], 'effect', group)
    
    scores_by_group[[i]][,A1:=a1]
    
    scores_by_group[[i]] <- scores_by_group[[i]][!is.na(a1)]
    scores_by_group[[i]][,c('a1','a2'):=NULL]
    
}

scores_by_group <- scores_by_group[lengths(scores_by_group) != 0]

# merge and sum up the predictor groups
merged <- Reduce(function(x,y){merge(x,y,by = c('SNP','A1'), all.x = T, all.y = T)}, scores_by_group)
merged$MultiPRS <- rowSums(merged[,3:ncol(merged)], na.rm = T)
merged <- merged[MultiPRS != 0]

cat(paste0('Final score contains ', nrow(merged), ' variants.\n'))

outfile <- gsub('(model)?\\.rds','multiprs.score.gz', opt$model_rds)
cat(paste0('Writing final score to ', outfile, '\n'))

fwrite(merged[,list(SNP, A1, MultiPRS)], outfile, sep='\t', quote=FALSE)