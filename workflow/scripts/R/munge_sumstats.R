
# A script that converts summary statistics to standard format

library(data.table)

# https://neurogenomics.github.io/MungeSumstats/articles/MungeSumstats.html

genome_version = '37'

sumstatsfile <- './METAANALYSIS_DIAGRAM_SE1_reformat.txt'

if (!require("BiocManager"))
  install.packages("BiocManager")

BiocManager::install("MungeSumstats")


if ( genome_version == '37'){
  # GRCh37
  BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
  BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")  
} else if (genome_version == '38') {
  # GRCh38
  BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh38")
  BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
} else {
  stop(paste0('Genome version "', genome_version, '" not supported'))
}

reformatted <- MungeSumstats::format_sumstats(path=sumstatsfile,ref_genome=paste0("GRCh", genome_version))

ss <- fread(reformatted, sep='\t') # this can be a huge table, and loading it this way might be a bad idea some times...

# renaming of non-standard columns
# for example: 
# setnames(ss, 'EAF_HAPMAP_CEU', 'FRQ')

fwrite(ss, file='METAANALYSIS_DIAGRAM_SE1_reformat.munged.txt.gz', sep = '\t')



