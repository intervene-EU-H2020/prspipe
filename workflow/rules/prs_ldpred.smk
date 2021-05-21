# Rules for preparing score and scale files for polygenic scoring using ldpred

rule ldpred_prep:
    # Implements the ldpred method
    # Note 1: the ldpred software requires Python dependencies, specified in the ldpred.yaml conda environment
    # Note 2: because the R script uses the ldpred software scripts (instead of the Python package), the script exports the module to the pythonpath
    input:
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        super_pop_keep=rules.create_ancestry.output['super_pop_keep']
    output:
        "{}/Score_files_for_polygenic/ldpred/{{study}}/1KGPhase3.w_hm3.{{study}}.{{ancestry}}.scale".format(config['Geno_1KG_dir'])
    log:
        "logs/prs_ldpred_{study}.{ancestry}.log"
    conda:
        "../envs/ldpred.yaml"
    shell:
        "export PYTHONPATH=$PYTHONPATH:$PWD/workflow/scripts/ldpred && ("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_LDPred/polygenic_score_file_creator_LDPred.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input.qc_stats} "
        "--plink {config[plink1_9]} "
        "--memory 20000 "
        "--n_cores 1 "
        "--ldpred \"python -m ldpred\" "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/ldpred/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input.super_pop_keep} "
        ") &> {log}"


rule all_ldpred_prep:
    # Run this rule to run the ldpred method
    # IMPORTANT NOTE you have to run the rules using the --conda-not-block-search-path-envvars option because of the pythonpath adjustment
    # TODO try to find a way to avoid this issue
    input: 
        expand("{}/Score_files_for_polygenic/ldpred/{{study.study_id}}/1KGPhase3.w_hm3.{{study.study_id}}.{{study.ancestry}}.scale".format(config['Geno_1KG_dir']), study=studies.itertuples())
