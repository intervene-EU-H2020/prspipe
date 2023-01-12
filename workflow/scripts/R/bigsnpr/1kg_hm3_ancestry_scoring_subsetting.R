library(dplyr)
library(bigsnpr)
library(data.table)

# extract the genetic variants that are provided by bigsnpr for ancestry projection that overlap the HapMap3-1kg variants

# paths, currently hard-coded.

stopifnot(file.exists("resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping.tsv.gz"))

ref_path <- "resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping.tsv.gz"

stopifnot(file.exists("resources/bigsnpr/ref_freqs.csv.gz"))
stopifnot(file.exists("resources/bigsnpr/projection.csv.gz"))

ref_freq_path <- "resources/bigsnpr/ref_freqs.csv.gz"
projection_path <- "resources/bigsnpr/projection.csv.gz"

ref <- fread(ref_path)

setnames(ref, c('a2', 'pos_hg19'), c('a0', 'pos'))
ref$beta <- 1

all_freq <- bigreadr::fread2(ref_freq_path)
projection <- bigreadr::fread2(projection_path)

matched <- snp_match(ref, all_freq[1:5])

ancestries <- c('AFR','AMR','EAS','EUR','SAS')

# matched$avg_freq <- rowMeans(all_freq[matched$`_NUM_ID_`,-(1:5)])
# avg_freq <- matched$avg_freq[1:10000]
# diagnostic plots
# for ( a in ancestries){
#     freq_ref <- matched[1:10000, paste0('MAF_', a)]
#     print(plot(freq_ref, avg_freq, col=ifelse(matched$beta[1:10000]==1, 'red', 'blue'), main = a, xlab = '1kg-hm3 freq', ylab='bigsnpr ref freq'))
# }

matched <- data.table(matched)

for ( a in ancestries){
    col <- paste0('MAF_', a)
    matched[,(col):=ifelse(beta == 1, get(col), 1-get(col))]
}

# diagnostic plots
# for ( a in ancestries){
#     freq_ref <- matched[1:10000, get(paste0('MAF_', a))]
#     print(plot(freq_ref, avg_freq, col=ifelse(matched$beta[1:10000]==1, 'red', 'blue'), main = a, xlab = '1kg-hm3 freq', ylab='bigsnpr ref freq'))
# }

correction <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)

results <- sapply(ancestries, function(a){
    
    col <- paste0('MAF_', a)
    res <- snp_ancestry_summary(
        freq = matched[[col]],
        info_freq_ref = all_freq[matched$`_NUM_ID_`, -(1:5)],
        projection = projection[matched$`_NUM_ID_`, -(1:5)],
        correction = correction
    )
    
})

write.table(results, 'resources/1kg/1KGPhase3.w_hm3.bigsnpr_ancestry_prop_estimates.tsv', row.names = T, col.names = T)

all_freq_1kg_hm3 <- all_freq[matched$`_NUM_ID_`,]
projection_1kg_hm3 <- projection[matched$`_NUM_ID_`,]

fwrite(all_freq_1kg_hm3, 'resources/bigsnpr/ref_freqs.1kg.w_hm3.tsv', sep='\t')
fwrite(projection_1kg_hm3, 'resources/bigsnpr/projection.1kg.w_hm3.tsv', sep='\t')

system('gzip resources/bigsnpr/ref_freqs.1kg.w_hm3.tsv')
system('gzip resources/bigsnpr/projection.1kg.w_hm3.tsv')