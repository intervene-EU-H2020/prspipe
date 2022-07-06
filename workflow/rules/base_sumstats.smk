# Rules for downloading and running QC for the base data from the GWAS catalog
# Automatically downloads summary statistics for studies specified in config file

rule download_sumstats:
    # downloads summary statistics from the GWAS catalog for the studies specified in the studies.tsv file
    # If a local path is specified, it tries to rename the columns and copies the data to a new location.
    output:
        temp("resources/sumstats/{study}.{ancestry}.raw")
    log:
        "logs/download_sumstats/{study}_{ancestry}.log"
    resources:
        mem_mb=20000
    shell:
        "("
        "python workflow/scripts/Python/gwas_catalog_sumstats.py "
        "--study-id {wildcards[study]} "
        "--out {output} "
        "--studies {config[studies]} "
        ") &> {log} "
        
        
rule map_sumstats:
    # forces mapping of sumstats to hg19 , and intersection with HapMap3 v3 (using pre-defined liftover-derived coordinates)
    input:
        sumstats_raw="resources/sumstats/{study}.{ancestry}.raw",
        mapping="resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping.tsv.gz"
    log:
        "logs/map_sumstats/{study}.{ancestry}.log"
    output:
        sumstats="resources/sumstats/{study}.{ancestry}.w_hm3.tsv.gz"
    params:
        out_prefix = lambda wc, output: output['sumstats'][0:-3]
    resources:
        mem_mb=20000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh --no-container-mount-home"
    shell:
        "("
        "Rscript workflow/scripts/R/setup/sumstat_harmonizer.R "
        "--sumstats {input[sumstats_raw]} "
        "--map {input[mapping]} "
        "--force TRUE "
        "--output {params[out_prefix]} "
        "--gz TRUE "
        "--tmpdir ./temp "
        ") &> {log}"
        
        
# the input file of the rule below changes depending on the force_sumstats_remapping parameter
remap_sumstats = bool(config['force_sumstats_remapping']) if 'force_sumstats_remapping' in config else False
if remap_sumstats:
    qc_sumstats_input = "resources/sumstats/{study}.{ancestry}.w_hm3.tsv.gz"
else:
    qc_sumstats_input = "resources/sumstats/{study}.{ancestry}.raw"


rule QC_sumstats:
    # run the cleaner script for summary statistics QC
    input:
        sumstats=qc_sumstats_input,
        other=rules.run_allele_freq_superpop.output
    output:
        "resources/sumstats/{study}.{ancestry}.cleaned.gz"
    log:
        "logs/QC_sumstats/{study}_{ancestry}.log"
    singularity:
        config['singularity']['all']
    resources:
        mem_mb=20000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh --no-container-mount-home"
    shell:
        "("
        "{config[Rscript]} {config[GenoPred_dir]}/Scripts/sumstat_cleaner/sumstat_cleaner.R "
        "--sumstats {input[sumstats]} "
        "--ss_freq_col FRQ "
        "--ref_plink_chr resources/1kg/1KGPhase3.w_hm3.chr "
        "--ref_freq_chr resources/1kg/freq_files/{wildcards[ancestry]}/1KGPhase3.w_hm3.{wildcards[ancestry]}.chr "
        "--output {output}"
        ") &> {log}"


rule all_QC:
    # download summary statistics and run all QC
    input: 
        expand("resources/sumstats/{study}.{ancestry}.cleaned.gz", zip, study=studies.study_id, ancestry=studies.ancestry)
