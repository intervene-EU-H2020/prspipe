
###############################################################
# nested sparse thresholding i.e. pT + clump for 1000 Genomes #
###############################################################

# plink2 version should replace plink1 version.

rule prs_scoring_ptclump_sparse:
    # 4.2.2 Sparse thresholding (not nested)
    # Here we will only use 8 p-value thresholds.
    # This section uses an R script called ‘polygenic_score_file_creator.R’
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True)
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0]
    output:
        touch('prs/pt_clump/{study}/ok'),
        scale=expand("prs/pt_clump/{{study}}/1KGPhase3.w_hm3.{{study}}.{ancestry}.scale", ancestry=config['1kg_superpop']),
        score="prs/pt_clump/{study}/1KGPhase3.w_hm3.{study}.score.gz",
        range_list="prs/pt_clump/{study}/1KGPhase3.w_hm3.{study}.range_list"
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_3.sqsh --no-container-mount-home"
    log:
        "logs/prs_scoring_ptclump_sparse/{study}.log"
    threads:
        1
    shell:
        "( "
        "{config[Rscript]} {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator/polygenic_score_file_creator_plink2.R "
        "--ref_plink_chr resources/1kg/1KGPhase3.w_hm3.chr "
        "--ref_keep resources/1kg/keep_files/{params[study_ancestry]}_samples.keep "
        "--sumstats {input[qc_stats]} "
        "--plink1 {config[plink1_9]} "
        "--plink2 {config[plink2]} "
        "--nested F "
        "--memory 3000 "
        "--output prs/pt_clump/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--tmpdir ./temp"
        ") &> {log} "


rule all_prs_scoring_ptclump_sparse:
    # runs rule above for all studies
    input:
        expand(rules.prs_scoring_ptclump_sparse.output, study=studies.study_id)
