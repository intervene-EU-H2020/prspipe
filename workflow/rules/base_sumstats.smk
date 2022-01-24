# Rules for downloading and running QC for the base data from the GWAS catalog
# Automatically downloads summary statistics for studies specified in config file

rule download_sumstats:
    # downloads summary statistics from the GWAS catalog
    # for the studies specified in the studies.tsv file
    # and puts them in a format suitable for the QC script
    output:
        "resources/sumstats/{study}.{ancestry}.raw",
    log:
        "logs/download_sumstats/{study}_{ancestry}.log"
    shell:
        "python workflow/scripts/Python/gwas_catalog_sumstats.py "
        "--study-id {wildcards[study]} "
        "--out {output} "
        "--studies {config[studies]} "
        "| tee -a {log}"


rule QC_sumstats:
    # run the cleaner script for summary statistics QC
    input:
        sumstats="resources/sumstats/{study}.{ancestry}.raw",
        other=rules.run_allele_freq_superpop.output
    output:
        "resources/sumstats/{study}.{ancestry}.cleaned.gz"
    log:
        "logs/QC_sumstats/{study}_{ancestry}.log"
    singularity:
        config['singularity']['all']
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
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
