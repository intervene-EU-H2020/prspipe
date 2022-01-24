
rule prs_scoring_lassosum:
    # Implements the lassosum method
    input:
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        hm3_gw=rules.extract_hm3_gw.output
    output:
        touch('prs/lassosum/{study}/ok')
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0]
    log:
        "logs/prs_scoring_lassosum/{study}.log"
    threads:
        1
    resources:
        mem_mb=5000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    shell:
        "("
        "{config[Rscript]} {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_lassosum/polygenic_score_file_creator_lassosum_plink2.R "
        "--ref_plink_gw resources/1kg/1KGPhase3.w_hm3.GW "
        "--ref_keep resources/1kg/keep_files/{params[study_ancestry]}_samples.keep "
        "--sumstats {input.qc_stats} "
        "--output prs/lassosum/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--plink2 {config[plink2]} "
        "--memory {resources[mem_mb]} "
        "--ref_pop_scale {input.super_pop_keep} "
        ") &> {log}"


rule all_lassosum_prep:
    # Run this rule to run the lassosum method
    input: 
       expand(rules.prs_scoring_lassosum.output, study=studies.study_id)
