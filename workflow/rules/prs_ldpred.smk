# Rules for preparing score and scale files for polygenic scoring using ldpred

# TODO need to set ld_pred in config using https://github.com/bvilhjal/ldpred

rule ldpred_prep:
    input: 
        "{}/{{study}}.{{ancestry}}.cleaned.gz".format(config['Base_sumstats_dir'])
    output:
        "{}/ldpred/1KGPhase3.w_hm3.{{ancestry}}.{{study}}".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/prs_ldpred_{study}.{ancestry}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_LDPred/polygenic_score_file_creator_LDPred.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input} "
        "--plink {config[plink1_9]} "
        "--memory 20000 "
        "--n_cores 1 "
        "--ldpred {config[ld_pred]} "
        "--output {output} "
        "--ref_pop_scale {config['Geno_1KG_dir']}/super_pop_keep.list "
        ") &> {log}"


rule all_ldpred_prep:
    input: 
        expand("{}/ldpred/1KGPhase3.w_hm3.{{study.ancestry}}.{{study.study_id}}".format(config['Base_sumstats_dir']), study=studies.itertuples())
