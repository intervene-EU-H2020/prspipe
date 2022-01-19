
# TODO: add dependencies bigsnpr, bigreadr, runonce


rule download_ld_reference_ldpred2:
    # download EUR (UKBB) LD reference
    output:
        ld_chr=expand('resources/LD_matrix/ldpred2/UKBB/precomputed/EUR/LD_chr{chr}.rds', chr=range(1,23)),
        map='resources/LD_matrix/ldpred2/UKBB/precomputed/EUR/map.rds'
    log:
        "logs/download_ld_reference_ldpred2.log"
    shell:
        "("
        "cd \"$(dirname {output[map]})\"; "
        "wget -O ldpred2_reference.zip https://ndownloader.figshare.com/articles/13034123/versions/3 ; "
        "unzip ldpred2_reference.zip && rm ldpred2_reference.zip ; "
        ") &> {log} "


rule run_ldpred2_precompld_1kg_deprecated:
    # Deprecated.
    # use plink2-version below instead!
    # only still here because the (also deprecated) ukbb.smk rules depend on it...
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        ld_ref=rules.download_ld_reference_ldpred2.output
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        is_binary=lambda wc: {'no':'FALSE', 'yes':'TRUE'}[studies.binary[studies.study_id == wc.study].iloc[0]]
    output:
        # score_grid=expand('{geno1kg}/Score_files_for_polygenic/LDPred2_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.{s1}_{s2}_{sparse}.SCORE', geno1kg=config['Geno_1KG_dir'], s1=suff1, s2=suff2, sparse=['sparse','nosparse']),
        # score_inf='{}/Score_files_for_polygenic/LDPred2_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.beta_inf.SCORE'.format(config['Geno_1KG_dir']),
        #log='{}/Score_files_for_polygenic/LDPred2_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.log'.format(config['Geno_1KG_dir']),
        scale=expand('{geno1kg}/Score_files_for_polygenic/LDPred2_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.{superpop}.scale', geno1kg=config['Geno_1KG_dir'], superpop=config['1kg_superpop'])
    log:
        "logs/run_ldpred2_precompld_1kg_deprecated_{study}.log"
    threads:
        16
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_LDPred2/polygenic_score_file_creator_LDPred2.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{params[study_ancestry]}_samples.keep "
        "--ldpred2_ref_dir resources/LD_matrix/ldpred2/UKBB/precomputed/EUR "
        "--sumstats {input[qc_stats]} "
        "--plink {config[plink1_9]} "
        "--memory 20000 "
        "--n_cores {threads} "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/LDPred2_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--ldpred2_ref_precomputed TRUE "
        "--binary {params[is_binary]} "
        ") &> {log} "


rule run_ldpred2_precompld_1kg_plink2:
    # note: this will break if GWAS ancestry is not EUR!
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        ld_ref=rules.download_ld_reference_ldpred2.output,
        extr_hm3_chr=expand(rules.extract_hm3.output, chr=range(1,23)),
        extr_hm3_gw=rules.extract_hm3_gw.output
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        is_binary=lambda wc: {'no':'FALSE', 'yes':'TRUE', 'No':'FALSE', 'Yes':'TRUE'}[studies.binary[studies.study_id == wc.study].iloc[0]]
    output:
        touch('prs/ldpred2/{study}/ok')
    log:
        "logs/run_ldpred2_precompld_1kg_{study}.log"
    threads:
        10
    resources:
        mem_mb=80000,
        time="14:00:00",
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    shell:
        "( "
        "export OPENBLAS_NUM_THREADS=1; "
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_LDPred2/polygenic_score_file_creator_LDPred2_LDPredRef_plink2.R "
        "--ref_plink resources/1kg/1KGPhase3.w_hm3.GW "
        "--ref_keep resources/1kg/keep_files/{params[study_ancestry]}_samples.keep "
        "--ldpred2_ref_dir resources/LD_matrix/ldpred2/UKBB/precomputed/EUR "
        "--sumstats {input[qc_stats]} "
        "--plink2 {config[plink2]} "
        "--memory {resources.mem_mb} "
        "--n_cores {threads} "
        "--output prs/ldpred2/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--binary {params[is_binary]} || "
        "rm prs/ldpred2/{wildcards[study]}/*bk prs/ldpred2/{wildcards[study]}/*rds "
        ") &> {log} "
        
        
rule all_run_ldpred2_1kg_precompld_1kg:
    input:
        expand(rules.run_ldpred2_precompld_1kg_plink2.output, study=studies.study_id)

