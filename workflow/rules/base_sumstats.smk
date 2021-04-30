# Rules for downloading and running QC for the base data from the GWAS catalog
# Automatically downloads summary statistics for studies specified in config file


rule download_sumstats:
    # downloads summary statistics from the GWAS catalog
    # for the specific study ids 
    output:
        "{}/{{study}}.raw".format(config['Base_sumstats_dir']),
    log:
        "logs/base_sumstats/download_{study}.log"
    shell:
        "("
        "python {config[Python_scr_dir]}/gwas_catalog_sumstats.py "
        "--study-id {wildcards[study]} "
        "--out {output}"
        ") &> {log}"


rule QC_sumstats:
    # run the cleaner script for summary statistics QC
    # TODO at the moment this only runs QC for EUR ancestry, but can generalize to allow user to input whatever they want
    input:
        "{}/{{study}}.raw".format(config['Base_sumstats_dir'])
    output:
        "{}/{{study}}.cleaned".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/qc_{study}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/sumstat_cleaner/sumstat_cleaner.R "
        "--sumstats {input} "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--ref_freq_chr {config[Geno_1KG_dir]}/freq_files/EUR/1KGPhase3.w_hm3.EUR.chr "
        "--gz F "
        "--output {output}"
        ") &> {log}"

# TODO the above rules implement step 4.1 of GenoPred - 
# the next steps from here are to create rules for the scripts
# implementing pre-processing steps for each polygenic score method

rule all_QC:
    # download summary statistics and run all QC
    input: 
        expand(rules.QC_sumstats.output, study=config['studies'])

