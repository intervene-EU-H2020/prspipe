
# rules to estimate ancestry proportions from summary statitics data using bigsnpr. 
# the scripts are not very refined, and use hard-coded paths.
# defining a conda-environment by its name has only become possible with more recent snakemake versions.
# the conda environment "bigsnpr" I used had the following R-packages installed:

"""

R version 4.2.2 (2022-10-31)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dplyr_1.0.10      bigsnpr_1.11.6    bigstatsr_1.5.12  data.table_1.14.4

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.9         pillar_1.8.1       compiler_4.2.2     bigsparser_0.6.1  
 [5] bigparallelr_0.3.2 iterators_1.0.14   rngtools_1.5.2     digest_0.6.30     
 [9] lifecycle_1.0.3    tibble_3.1.8       gtable_0.3.1       lattice_0.20-45   
[13] pkgconfig_2.0.3    rlang_1.0.6        doRNG_1.8.2        Matrix_1.5-3      
[17] foreach_1.5.2      cli_3.4.1          flock_0.7          parallel_4.2.2    
[21] bigassertr_0.1.5   generics_0.1.3     vctrs_0.5.1        grid_4.2.2        
[25] tidyselect_1.2.0   cowplot_1.1.1      glue_1.6.2         R6_2.5.1          
[29] fansi_1.0.3        ggplot2_3.4.0      magrittr_2.0.3     scales_1.2.1      
[33] codetools_0.2-18   colorspace_2.0-3   utf8_1.2.2         munsell_0.5.0     
[37] doParallel_1.0.17

"""


rule bigsnpr_download_allele_freq_projections:
    output:
        ref_freqs = "resources/bigsnpr/ref_freqs.csv.gz",
        projection = "resources/bigsnpr/projection.csv.gz"
    shell:
        "wget -O resources/bigsnpr/ref_freqs.csv.gz https://figshare.com/ndownloader/files/31620968 && "
        "wget -O resources/bigsnpr/projection.csv.gz https://figshare.com/ndownloader/files/31620953"


rule biobank_allele_freq_ancestry:
    # needs plink2
    input:
        geno=expand('custom_input/{bbid}/genotypes/chr{chrom}.{ext}', chrom=range(1,23), ext=['bed','bim','fam'], allow_missing=True),
        keepfiles=expand("results/{bbid}/Ancestry_identifier/outlier_detection/AllAncestry.QC.{superpop}.keep", superpop=config['1kg_superpop'], allow_missing=True)
    output:
        acount = expand('results/{bbid}/afreq/AllAncestry.QC.{superpop}_chr{chrom}.acount', superpop=config['1kg_superpop'], chrom=range(1,23), allow_missing=True),
        vmiss = expand('results/{bbid}/afreq/AllAncestry.QC.{superpop}_chr{chrom}.vmiss', superpop=config['1kg_superpop'], chrom=range(1,23), allow_missing=True)
    shell:
        "export BIOBANK={wildcards[bbid]} && "
        "bash workflow/scripts/bash/export_allele_frequencies_missingness.sh"
        
        
rule bigsnpr_hm3_1kg_ref_ancestry_proportions:
    # needs bigsnpr
    input:
        "resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping.tsv.gz",
        "resources/bigsnpr/ref_freqs.csv.gz",
        "resources/bigsnpr/projection.csv.gz"
    output:
        ancestry_prop_estimates="resources/1kg/1KGPhase3.w_hm3.bigsnpr_ancestry_prop_estimates.tsv",
        ref_freqs = "resources/bigsnpr/ref_freqs.1kg.w_hm3.tsv.gz",
        projection = "resources/bigsnpr/projection.1kg.w_hm3.tsv"
    conda:
        "bigsnpr"
    log:
        "logs/bigsnpr_hm3_1kg_ref_ancestry_proportions.log"
    shell:
        "("
        "Rscript workflow/scripts/R/bigsnpr/1kg_hm3_ancestry_scoring_subsetting.R "
        ") &> {log}"
        
        
rule bigsnpr_biobank_ancestry_proportions:
    # needs bigsnpr
    input:
        ref_mapping="resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping.tsv.gz",
        ref_freqs="resources/bigsnpr/ref_freqs.1kg.w_hm3.tsv.gz",
        projection="resources/bigsnpr/projection.1kg.w_hm3.tsv.gz",
        acount=rules.biobank_allele_freq_ancestry.output['acount'],
        vmiss=rules.biobank_allele_freq_ancestry.output['vmiss']
    output:
        bigsnpr_ancestry_prop_results="results/{bbid}/afreq/AllAncestry.QC.{superpop}.bigsnpr_ancestry_prop_results.tsv",
        diagnostic_plot="results/{bbid}/afreq/AllAncestry.QC.{superpop}.freq_scatter.png"
    params:
        prefix=lambda wc: f"results/{wc['bbid']}/afreq/AllAncestry.QC.{wc['superpop']}_chr"
    conda:
        "bigsnpr"
    log:
        "logs/bigsnpr_biobank_ancestry_proportions/{bbid}/{superpop}.log"
    shell:
        "("
        "Rscript workflow/scripts/R/bigsnpr/ancestry_from_afreq.R {params[prefix]}"
        ") &> {log}"
        
        
        
rule all_bigsnpr_biobank_ancestry_proportions:
    # rule that requests all the outputs for a specific biobank
    input:
        expand(rules.bigsnpr_biobank_ancestry_proportions.output, superpop=config['1kg_superpop'], allow_missing=True)
    output:
        temp(touch("results/{bbid}/afreq/biobank_ancestry_proportions.ok"))
        
        
        
###########################################################
# check the overlap of scoring files and target genotypes #
###########################################################


rule score_to_target_geno_overlap_prep:
    # create a temporary file to load in the rule below
    input:
        expand('prs/{method}/{study}/1KGPhase3.w_hm3.{study}.score.gz', method=config['prs_methods'], study=studies.study_id)
    output:
        temp('temp/scores.txt')
    run:
        with open(output[0], 'w') as outfile:
            for f in input:
                outfile.write(f'{f}\n')

localrules:
    score_to_target_geno_overlap_prep
    
    
rule score_to_target_geno_overlap_create_rsid_rds:
    # export all the rsids into a single file
    input:
        'temp/scores.txt'
    output:
        'temp/prs_rsid.rds'
    log:
        'logs/score_to_target_geno_create_rsid_rds.log'
    threads:
        1
    resources:
        mem_mb=8000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home",
        time="1:00:00"
    shell:
        "("
        "Rscript workflow/scripts/R/get_rsid_from_score_files.R {input}"
        ") &> {log} "
        
        

rule score_to_target_geno_overlap_calculate_overlaps:
    input:
        rsid_rds='temp/prs_rsid.rds',
        acount=rules.biobank_allele_freq_ancestry.output['acount'],
        vmiss=rules.biobank_allele_freq_ancestry.output['vmiss']
    output:
        bigsnpr_ancestry_prop_results="results/{bbid}/afreq/AllAncestry.QC.{superpop}.score_overlaps.tsv.gz"
    params:
        prefix=lambda wc: f"results/{wc['bbid']}/afreq/AllAncestry.QC.{wc['superpop']}_chr"
    conda:
        "bigsnpr"
    log:
        "logs/score_to_target_geno_overlap_calculate_overlaps/{bbid}/{superpop}.log"
    shell:
        "("
        "Rscript workflow/scripts/R/export_score_file_target_overlap.R {params[prefix]}"
        ") &> {log}"