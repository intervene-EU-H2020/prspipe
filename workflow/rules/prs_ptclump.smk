

rule nested_sparse_thresholding_1kg:
    # 4.2.1 https://opain.github.io/GenoPred/Pipeline_prep.html
    # Here I prepare reference files for typical polygenic scores derived using the p-value thresholding and LD-based clumping procedure.
    # Here we will only use 8 p-value thresholds.
    # This section uses an R script called ‘polygenic_score_file_creator.R’
    # https://github.com/opain/GenoPred/blob/master/Scripts/polygenic_score_file_creator/polygenic_score_file_creator.R
    input:
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        super_pop_keep=rules.create_ancestry.output['super_pop_keep']
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0]
    output:
        scale=expand("{geno1kdir}/Score_files_for_polygenic/pt_clump/{{study}}/1KGPhase3.w_hm3.{{study}}.{ancestry}.scale", geno1kdir=[config['Geno_1KG_dir']], ancestry=config['1kg_superpop']),
        score=expand("{geno1kdir}/Score_files_for_polygenic/pt_clump/{{study}}/1KGPhase3.w_hm3.{{study}}.chr{chr_id}.score", geno1kdir=[config['Geno_1KG_dir']], chr_id=range(1,23)),
        range_values=expand("{geno1kdir}/Score_files_for_polygenic/pt_clump/{{study}}/1KGPhase3.w_hm3.{{study}}.chr{chr_id}.range_values", geno1kdir=[config['Geno_1KG_dir']], chr_id = range(1,23))
    log:
        "logs/nested_sparse_thresholding_1kg/{study}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator/polygenic_score_file_creator.R "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{params[study_ancestry]}_samples.keep "
        "--sumstats {input.qc_stats} "
        "--plink {config[plink1_9]} "
        "--memory 3000 "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input.super_pop_keep} "
        ") &> {log} "


rule all_nested_sparse_thresholding_1kg:
    # runs rule above for all studies
    input:
        expand(rules.nested_sparse_thresholding_1kg.output, study=studies.study_id)


rule sparse_thresholding_1kg:
    # 4.2.2 Sparse thresholding (not nested)
    # Here we will only use 8 p-value thresholds.
    # This section uses an R script called ‘polygenic_score_file_creator.R’
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True)
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0]
    output:
        scale=expand("{geno1kdir}/Score_files_for_polygenic/pt_clump_nonnested/{{study}}/1KGPhase3.w_hm3.{{study}}.{ancestry}.scale", geno1kdir=[config['Geno_1KG_dir']], ancestry=config['1kg_superpop']),
        score=expand("{geno1kdir}/Score_files_for_polygenic/pt_clump_nonnested/{{study}}/1KGPhase3.w_hm3.{{study}}.chr{chr_id}.score", geno1kdir=[config['Geno_1KG_dir']], chr_id=range(1,23)),
        range_values=expand("{geno1kdir}/Score_files_for_polygenic/pt_clump_nonnested/{{study}}/1KGPhase3.w_hm3.{{study}}.chr{chr_id}.range_values", geno1kdir=[config['Geno_1KG_dir']], chr_id = range(1,23))
    log:
        "logs/sparse_thresholding_1kg/{study}.log"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator/polygenic_score_file_creator.R "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{params[study_ancestry]}_samples.keep "
        "--sumstats {input[qc_stats]} "
        "--plink {config[plink1_9]} "
        "--nested F "
        "--memory 3000 "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump_nonnested/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input[super_pop_keep]} "
        ") &> {log} "


rule all_sparse_thresholding_1kg:
    # runs rule above for all studies
    input:
        expand(rules.sparse_thresholding_1kg.output, study=studies.study_id)


rule dense_thresholding_1kg:
    # 4.2.3 Dense thresholding
    # Here we will only use dense p-value thresholds, mimicking the PRSice approach.
    # This section uses an R script called ‘polygenic_score_file_creator.R’
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True)
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0]
    output:
        scale=expand("{geno1kdir}/Score_files_for_polygenic/pt_clump_dense/{{study}}/1KGPhase3.w_hm3.{{study}}.{ancestry}.scale", geno1kdir=[config['Geno_1KG_dir']], ancestry=config['1kg_superpop']),
        # score=expand("{geno1kdir}/Score_files_for_polygenic/pt_clump_dense/{{study}}/1KGPhase3.w_hm3.{{study}}.chr{chr_id}.score", geno1kdir=[config['Geno_1KG_dir']], chr_id=range(1,23)),
        # range_values=expand("{geno1kdir}/Score_files_for_polygenic/pt_clump_dense/{{study}}/1KGPhase3.w_hm3.{{study}}.chr{chr_id}.range_values", geno1kdir=[config['Geno_1KG_dir']], chr_id = range(1,23))
    log:
        "logs/dense_thresholding_1kg/{study}.log"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_dense/polygenic_score_file_creator_dense.R "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{params[study_ancestry]}_samples.keep "
        "--sumstats {input[qc_stats]} "
        "--plink {config[plink1_9]} "
        "--prsice_path {config[prsice_path]} "
        "--memory 3000 "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump_dense/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input[super_pop_keep]} "
        ") &> {log} "


rule all_dense_thresholding_1kg:
    # runs rule above for all studies
    input:
        expand(rules.dense_thresholding_1kg.output, study=studies.study_id)
