# Rules for preparing score and scale files for polygenic scoring using lassosum

rule lassosum_prep:
    # Implements the lassosum method
    input: 
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        super_pop_keep=rules.create_ancestry.output['super_pop_keep']
    output:
        "{}/Score_files_for_polygenic/lassosum/{{study}}/1KGPhase3.w_hm3.{{study}}.{{ancestry}}.scale".format(config['Geno_1KG_dir'])
    log:
        "logs/prs_lassosum_{study}.{ancestry}.log"
    threads:
        1
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_lassosum/polygenic_score_file_creator_lassosum.R "
        "--ref_plink_gw {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input.qc_stats} "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/lassosum/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--plink {config[plink1_9]} "
        "--ref_pop_scale {input.super_pop_keep} "
        ") &> {log}"


rule all_lassosum_prep:
    # Run this rule to run the lassosum method
    input: 
        expand("{}/Score_files_for_polygenic/lassosum/{{study.study_id}}/1KGPhase3.w_hm3.{{study.study_id}}.{{study.ancestry}}.scale".format(config['Geno_1KG_dir']), study=studies.itertuples())
