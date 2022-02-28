#!/usr/bin/enb Rscript

# script derived from the GenoPred pipeline
# https://opain.github.io/GenoPred/
# original author: Oliver Pain
# adapted by Remo Monti

# Here we download information on which population each individual in the 1000 Genomes reference is from. Individuals are grouped into ‘populations’ which are typically country specific, and ‘super populations’ which include a collection of ancetrally similar countries. Individuals are grouped into these populations if the last few generations of their family are all from one region. We need this population data so we can select individuals in the 1000 Genomes data to match those in our target samples. This is important for providing accurate information on LD structure and minor allele frequencies.

# Read in environmental variables
source('workflow/scripts/R/source_config.R')

pop_data<-read.table(paste0('resources/1kg/integrated_call_samples_v3.20130502.ALL.panel'), header=T, stringsAsFactors=F)

for(i in unique(pop_data$pop)){
  write.table(cbind(pop_data$sample[pop_data$pop == i],pop_data$sample[pop_data$pop == i]), paste0('resources/1kg/keep_files/',i,'_samples.keep'), col.names=F, row.names=F, quote=F)
}

for(i in unique(pop_data$super_pop)){
  write.table(cbind(pop_data$sample[pop_data$super_pop == i],pop_data$sample[pop_data$super_pop == i]), paste0('resources/1kg/keep_files/',i,'_samples.keep'), col.names=F, row.names=F, quote=F)
}

# Create a file listing all ancestry groups
write.table(unique(pop_data$super_pop), paste0('resources/1kg/super_pop.list'), col.names=F, row.names=F, quote=F)
write.table(unique(pop_data$pop), paste0('resources/1kg/pop.list'), col.names=F, row.names=F, quote=F)

# Create a file listing the code of each population and the location of the keep file
pop_keep_loc<-data.frame(pop=unique(pop_data$pop))
pop_keep_loc$keep<-paste0('resources/1kg/keep_files/',pop_keep_loc$pop,'_samples.keep')

super_pop_keep_loc<-data.frame(pop=unique(pop_data$super_pop))
super_pop_keep_loc$keep<-paste0('resources/1kg/keep_files/',super_pop_keep_loc$pop,'_samples.keep')

write.table(super_pop_keep_loc, paste0('resources/1kg/super_pop_keep.list'), col.names=F, row.names=F, quote=F)
write.table(pop_keep_loc, paste0('resources/1kg/pop_keep.list'), col.names=F, row.names=F, quote=F)
write.table(rbind(super_pop_keep_loc,pop_keep_loc), paste0('resources/1kg/super_pop_and_pop_keep.list'), col.names=F, row.names=F, quote=F)
