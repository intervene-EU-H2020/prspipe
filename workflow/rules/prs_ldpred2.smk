
# TODO: add dependencies bigsnpr, bigreadr, runonce

rule generate_ld_referece_ldpred2:
    # 4.8.1 Create LD reference for LDPred2, part 1
    # uses a lot of memory > 16G
    # Mostly single core seems to be used, adding cores probably won't speed it up much
    # TODO: move parameters to config (?)
    input:
        keep=rules.create_ancestry.output['super_pop_keep'],
        genetic_map = rules.download_genetic_map.output,
        bed_bim_fam = rules.extract_hm3.output
    output:
        rds=expand('resources/LD_matrix/ldpred2/1000G/fromscratch/{popul}/LD_chr{chr}.rds', chr=range(1,23), allow_missing=True),
        map='resources/LD_matrix/ldpred2/1000G/fromscratch/{popul}/map.rds',
        sd='resources/LD_matrix/ldpred2/1000G/fromscratch/{popul}/sd.rds'
    params:
        cores=6
    log:
        "logs/generate_ld_referece_ldpred2_{popul}.log"
    shell:
        "("
        "Rscript workflow/scripts/R/prs_ldpred2/generate_ld_reference.R {wildcards[popul]} {params[cores]} "
        ") &> {log} "
        
        
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


#rule generate_sd_file_ldpred2:
#    # 4.8.1 Create LD reference for LDPred2, part 2
#    # RM: not needed. 
#    input:
#        map=rules.generate_ld_referece_ldpred2.output.map
#    output:
#        sd='{}/LD_matrix/LDPred2/{{popul}}/sd.rds'.format(config['Geno_1KG_dir']),
#        map='{}/LD_matrix/LDPred2/{{popul}}/map.rds'.format(config['Geno_1KG_dir'])
#    log:
#        "logs/generate_sd_file_ldpred2_{popul}.log"
#    shell:
#        "("
#        "Rscript workflow/scripts/R/prs_ldpred2/generate_sd_file.R {input[map]} " 
#        ") &> {log}"


# TODO: make this work! (see "change outputs" below)
# suff1 = ["0.00018", "0.00032", "0.00056", "0.0018", "0.001", "0.0032", "0.0056", "0.018", "0.01", "0.032", "0.056", "0.18", "0.1", "0.32", "0.56", "1.8e.05", "1", "1e.04", "1e.05", "3.2e.05", "5.6e.05"]
# suff2 = ["0.0012","6e.04","8e.04"]

rule run_ldpred2_precompld_1kg:
    # note: this will break if GWAS ancestry is not EUR!
    # TODO: move parameters to config
    # TODO: check if phenotype is binary!
    # TODO: currently this assumes we only want predictions for the same ancestry as the GWAS was performed in
    # TODO: change outputs
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
        log='{}/Score_files_for_polygenic/LDPred2_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.log'.format(config['Geno_1KG_dir']),
        scale=expand('{geno1kg}/Score_files_for_polygenic/LDPred2_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.{superpop}.scale', geno1kg=config['Geno_1KG_dir'], superpop=config['1kg_superpop'])
        # touch('run_ldpred2_precompld_1kg_{study}.ok')
    log:
        "logs/run_ldpred2_precompld_1kg_{study}.log"
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
        
        
rule run_ldpred2_1kg:
    # TODO: move parameters to config
    # TODO: check if phenotype is binary!
    # TODO: currently this assumes we only want predictions for the same ancestry as the GWAS was performed in
    # TODO: change outputs
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        ld_ref=lambda wc: expand(rules.generate_ld_referece_ldpred2.output, popul=studies.ancestry[studies.study_id == wc.study], allow_missing=True)
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        is_binary=lambda wc: {'no':'FALSE', 'yes':'TRUE'}[studies.binary[studies.study_id == wc.study].iloc[0]]
    output:
        # score_grid=expand('{geno1kg}/Score_files_for_polygenic/LDPred2/{{study}}/1KGPhase3.w_hm3.{{study}}.{s1}_{s2}_{sparse}.SCORE', geno1kg=config['Geno_1KG_dir'], s1=suff1, s2=suff2, sparse=['sparse','nosparse']),
        # score_inf='{}/Score_files_for_polygenic/LDPred2/{{study}}/1KGPhase3.w_hm3.{{study}}.beta_inf.SCORE'.format(config['Geno_1KG_dir']),
        log='{}/Score_files_for_polygenic/LDPred2/{{study}}/1KGPhase3.w_hm3.{{study}}.log'.format(config['Geno_1KG_dir']),
        scale=expand('{geno1kg}/Score_files_for_polygenic/LDPred2/{{study}}/1KGPhase3.w_hm3.{{study}}.{superpop}.scale', geno1kg=config['Geno_1KG_dir'], superpop=config['1kg_superpop'])
    log:
        "logs/run_ldpred2_1kg_{study}.log"
    threads:
        16
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_LDPred2/polygenic_score_file_creator_LDPred2.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{params[study_ancestry]}_samples.keep "
        "--ldpred2_ref_dir resources/LD_matrix/ldpred2/1000G/fromscratch/{params[study_ancestry]}/ "
        "--sumstats {input[qc_stats]} "
        "--plink {config[plink1_9]} "
        "--memory 20000 "
        "--n_cores {threads} "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/LDPred2/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--ldpred2_ref_precomputed FALSE "
        "--binary {params[is_binary]} "
        ") &> {log} "
        
        
rule all_run_ldpred2_1kg:
    # runs rules above
    input:
        expand(rules.run_ldpred2_1kg.output, study=studies.study_id),
        expand(rules.run_ldpred2_precompld_1kg.output, study=studies.study_id)
