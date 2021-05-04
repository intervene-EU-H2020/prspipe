# Rules for downloading and running QC for the base data from the GWAS catalog
# Automatically downloads summary statistics for studies specified in config file


rule download_sumstats:
    # downloads summary statistics from the GWAS catalog
    # for the studies specified in the studies.tsv file
    # and puts them in a format suitable for the QC script
    output:
        "{}/{{study}}.{{ancestry}}.raw".format(config['Base_sumstats_dir']),
    log:
        "logs/base_sumstats/download_{study}.{ancestry}.log"
    conda:
        "envs/sumstats.yaml"
    shell:
        "python {config[Python_scr_dir]}/gwas_catalog_sumstats.py "
        "--study-id {wildcards[study]} "
        "--out {output} "
        "| tee -a {log}"


rule QC_sumstats:
    # run the cleaner script for summary statistics QC
    input:
        "{}/{{study}}.{{ancestry}}.raw".format(config['Base_sumstats_dir'])
    output:
        "{}/{{study}}.{{ancestry}}.cleaned.gz".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/qc_{study}.{ancestry}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/sumstat_cleaner/sumstat_cleaner.R "
        "--sumstats {input} "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--ref_freq_chr {config[Geno_1KG_dir]}/freq_files/{wildcards[ancestry]}/1KGPhase3.w_hm3.{wildcards[ancestry]}.chr "
        "--output {output}"
        ") &> {log}"


rule all_QC:
    # download summary statistics and run all QC
    input: 
        expand("{}/{{study.study_id}}.{{study.ancestry}}.cleaned.gz".format(config['Base_sumstats_dir']), study=studies.itertuples())