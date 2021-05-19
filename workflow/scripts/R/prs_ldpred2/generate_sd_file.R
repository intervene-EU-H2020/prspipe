
library(bigsnpr)
library(bigreadr)

# GenoPred 4.8.1, part 2
# Original author: Oliver Pain
# adapted by: Remo Monti

# this script is actually not used!

# "This was something else I considered but didn't finish. The second code block is me preparing the input files based on the reference data that is shared by developers of LDPred2. It doesn't overwrite the output from the previous code block, rather it is creating an alternative set of reference files. However, I never finished this component." - O.P.

# read configuration
source('workflow/scripts/R/source_config.R')

args = commandArgs(trailingOnly=TRUE)

# the input map file
infile = args[1]

out_map = paste0(dirname(infile), '/map.rds')
out_sd = paste0(dirname(infile), '/sd.rds')

# map<-readRDS('/users/k1806347/brc_scratch/Data/LDPred2/map.rds')
map <- readRDS(infile)
sd_ldref <- sqrt(2 * map$af_UKBB * (1 - map$af_UKBB))
# saveRDS(sd_ldref, '/users/k1806347/brc_scratch/Data/LDPred2/sd.rds')
print(head(sd_ldref))
saveRDS(sd_ldref, out_sd)

# Modify header in map
map$genetic.dist<-NA
print(head(map))
map<-map[c('chr','rsid','genetic.dist','pos','a1','a0','ld')]
names(map)<-c('chromosome','marker.ID','genetic.dist','physical.pos','allele1','allele2','ld')

saveRDS(map, out_map)