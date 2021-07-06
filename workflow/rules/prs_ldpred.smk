# Rules for preparing score and scale files for polygenic scoring using ldpred

rule ldpred_prep:
    # Implements the ldpred method (1000G reference)
    # Note 1: the ldpred software requires Python dependencies, specified in the ldpred.yaml conda environment
    input:
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        super_pop_keep=rules.create_ancestry.output['super_pop_keep']
    output:
        "{}/Score_files_for_polygenic/ldpred/{{study}}/1KGPhase3.w_hm3.{{study}}.{{ancestry}}.scale".format(config['Geno_1KG_dir'])
    log:
        "logs/prs_ldpred_{study}.{ancestry}.log"
    conda:
        "../envs/ldpred.yaml"
    threads:
        10
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_LDPred/polygenic_score_file_creator_LDPred.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input.qc_stats} "
        "--plink {config[plink1_9]} "
        "--memory 20000 "
        "--n_cores {threads} "
        "--ldpred ldpred "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/ldpred/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input.super_pop_keep} "
        ") &> {log}"


rule all_ldpred_prep:
    # Run this rule to run the ldpred method
    input: 
        expand("{}/Score_files_for_polygenic/ldpred/{{study.study_id}}/1KGPhase3.w_hm3.{{study.study_id}}.{{study.ancestry}}.scale".format(config['Geno_1KG_dir']), study=studies.itertuples())
