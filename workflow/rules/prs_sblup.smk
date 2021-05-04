# Rules for preparing score and scale files for polygenic scoring using sblup

# TODO need to set ldsc and ldsc_ref in config

rule sblup_prep:
    input: 
        "{}/{{study}}.{{ancestry}}.cleaned.gz".format(config['Base_sumstats_dir'])
    output:
        "{}/sblup/1KGPhase3.w_hm3.{{ancestry}}.{{study}}".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/prs_sblup_{study}.{ancestry}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_SBLUP/polygenic_score_file_creator_SBLUP.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input} "
        "--plink {config[plink1_9]} "
        "--gcta {config[gcta]} "
        "--munge_sumstats {config[GenoPred_dir]}/ldsc "
        "--ldsc {config[ldsc]} "
        "--ldsc_ref {config[ldsc_ref]} "
        "--hm3_snplist {config[HapMap3_snplist_dir]}/w_hm3.snplist "
        "--memory 50000 "
        "--n_cores 6 "
        "--output {output} "
        "--ref_pop_scale {config['Geno_1KG_dir']}/super_pop_keep.list "
        ") &> {log}"


rule all_sblup_prep:
    input: 
        expand("{}/sblup/1KGPhase3.w_hm3.{{study.ancestry}}.{{study.study_id}}".format(config['Base_sumstats_dir']), study=studies.itertuples())
