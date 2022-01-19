
# https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html
 
bbids = target_list['name'].values

############################
# Harmonization to HapMap3 #
############################

# See rules in genotype_harmonizatin.smk (harmonize_target_genotypes)


####################
# Ancestry scoring #
####################

# needs: 1000G genotypes, keep files, and integrated_call_samples_v3.20130502.ALL.panel_small

rule ancestry_scoring_ext:
    # ancestry scoring of custom genotypes
    # recommended 15G of RAM
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        harmonised_geno=expand("custom_input/{bbid}/genotypes/chr{chr}.bed", chr=range(1,23), allow_missing=True)
    output:
        eigenvec_sp=expand('results/{{bbid}}/Ancestry_idenitfier/AllAncestry.{superpop}.eigenvec', superpop=config['1kg_superpop']),
        keep_sp=expand('results/{{bbid}}/Ancestry_idenitfier/AllAncestry.{superpop}.keep', superpop=config['1kg_superpop']),
        scale_sp=expand('results/{{bbid}}/Ancestry_idenitfier/AllAncestry.{superpop}.scale', superpop=config['1kg_superpop']),
        eigenvec='results/{bbid}/Ancestry_idenitfier/AllAncestry.eigenvec',
        eigenvec_var='results/{bbid}/Ancestry_idenitfier/AllAncestry.eigenvec.var',
        log='results/{bbid}/Ancestry_idenitfier/AllAncestry.log',
        model_pred_keep=expand('results/{{bbid}}/Ancestry_idenitfier/AllAncestry.model_pred.{superpop}.keep', superpop=config['1kg_superpop']),
        model_pred='results/{bbid}/Ancestry_idenitfier/AllAncestry.model_pred',
        pop_model='results/{bbid}/Ancestry_idenitfier/AllAncestry.pop_model.rds',
        model_details='results/{bbid}/Ancestry_idenitfier/AllAncestry.pop_model_prediction_details.txt',
        scale='results/{bbid}/Ancestry_idenitfier/AllAncestry.scale'
    log:
        "logs/ancestry_scoring_ext/{bbid}.log"
    singularity:
        config['singularity']['all']
    resources:
        mem_mb=32000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home",
        time="03:00:00"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Ancestry_identifier/Ancestry_identifier.R "
        "--target_plink_chr custom_input/{wildcards[bbid]}/genotypes/chr "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--n_pcs 100 "
        "--plink {config[plink1_9]} "
        "--plink2 {config[plink2]} "
        "--output results/{wildcards[bbid]}/Ancestry_idenitfier/AllAncestry "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--pop_data {config[Geno_1KG_dir]}/integrated_call_samples_v3.20130502.ALL.panel_small "
        ") &> {log} "
        
       
rule calculate_maf_ancestry_ext:
    # calculate allele frequencies for different superpopulations
    # NA-values are dropped (TODO: make sure this doesn't lead to bugs downstream)
    input:
        keep_sp='results/{bbid}/Ancestry_idenitfier/AllAncestry.{superpop}.keep',
        harmonised_geno="custom_input/genotypes/{bbid}/chr{chr}.bed"
    output:
        afreq='results/{bbid}/Ancestry_idenitfier/AllAncestry.{superpop}.chr{chr}.frq'
    params:
        geno_prefix=lambda wc, input: input['harmonised_geno'][:-4],
        out_prefix=lambda wc, output: output['afreq'][:-4]
    log:
        'logs/calculate_maf_ancestry_ext/{bbid}/{superpop}_{chr}.log'
    singularity:
        config['singularity']['all']
    shell:
        '('
        '{config[plink1_9]} '
        '--bfile {params[geno_prefix]} '
        '--freq '
        '--keep {input[keep_sp]} '
        '--out {params[out_prefix]}_tmp ; '
        'awk \'$5 != "NA"\' {params[out_prefix]}_tmp.frq > {params[out_prefix]}.frq && rm {params[out_prefix]}_tmp.frq '
        ') &> {log} '
        # 'gzip {params[out_prefix]}.frq; ' 


rule all_calculate_maf_ancestry_ext:
    # runs rule above for all ancestries and chromosomes
    # Note: this fails if some ancestries are not represented, which is likely. Use "run.sh --keep-going ..." to force execution for all represented ancestries
    input:
        expand(rules.calculate_maf_ancestry_ext.output, chr=range(1,23), superpop=config['1kg_superpop'], bbid=bbids)
        

#################
# 10K reference #
#################

# we skip the generation of a 10K reference, see ukbb.smk if we wished to implement this for external biobanks as well.

        
#############################>>
# START: Polygenic scoring  #>>
#############################>>


# TODO:
# update this and make it use plink2 rules

# rule validate_setup_ext:
#     # requests all necessary outputs for the rules below.
#     # TODO: replace these with the plink2 counterparts
#     input:
#         # implemented by Sophie:
#         rules.all_dbslmm_prep.input,
#         rules.all_lassosum_prep.input,
#         rules.all_ldpred_prep.input,
#         rules.all_sblup_prep.input,
#         # implemented by Remo:
#         rules.all_run_prscs_precompld_1kg_refukbb.input,
#         # rules.all_run_ldpred2_1kg_precompld_1kg_.input
#         rules.all_sparse_thresholding_1kg.input,
#         rules.all_run_sbayesr_precompld_1kg_refukbb_robust.input


wildcard_constraints:
    score_id="[A-Za-z\\.0-9_\\-]+"

    
rule run_scaled_polygenic_scorer:
    # general purpose rule to run scaled_polygenic_scorer
    input:
        score='prs/{method}/{study}/{score_id}.score.gz',
        scale='prs/{method}/{study}/{score_id}.{superpop}.scale',
        ref_freq_chr=expand('resources/1kg/freq_files/{{superpop}}/1KGPhase3.w_hm3.{{superpop}}.chr{chr}.frq', chr=range(1,23)),
        target_plink=expand('custom_input/{{bbid}}/genotypes/chr{chr}.{suffix}', chr=range(1,23), suffix=['bed','bim','fam']),
        target_keep="results/{bbid}/Ancestry_idenitfier/AllAncestry.model_pred.{superpop}.keep"
    output:
        ok=touch('results/{bbid}/prs/{method}/{study}/{superpop}/{score_id}.{superpop}.done'),
        profiles='results/{bbid}/prs/{method}/{study}/{superpop}/{score_id}.{superpop}.profiles'
    params:
        geno_prefix=lambda wc, input: input['target_plink'][0][:-5],
        freq_prefix=lambda wc, input: input['ref_freq_chr'][0][:-5],
        out_prefix=lambda wc, output: output['ok'].replace('.done','')
    resources:
        mem_mb=8000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home",
        time="03:00:00"
    shell:
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer/Scaled_polygenic_scorer_plink2.R "
        "--target_plink_chr {params.geno_prefix} "
        "--target_keep {input.target_keep} "
        "--ref_score {input.score} "
        "--ref_scale {input.scale} "
        "--ref_freq_chr {params.freq_prefix} "
        "--plink2 {config[plink2]} "
        "--pheno_name {wildcards.study} "
        "--output {params[out_prefix]} "

        
rule all_run_scaled_polygenic_scorer:
    input:
        expand(rules.run_scaled_polygenic_scorer.output,
               bbid=target_list.name,
               study=studies.study_id.values,
               superpop=['EUR'],
               score_id=['1KGPhase3.w_hm3.{}'.format(s) for s in studies.study_id.values],
               method=['pt_clump','dbslmm'])



##########################
# Pruning & Thresholding #
##########################

# Deprecated rules. Use run_scaled_polygenic_scorer with method = "pt_clump" instead.

rule nested_sparse_thresholding_score_ext_ref1kg:
    # P+T clump - nested sparse
    # TODO: other ancestries 
    # TODO: move input files (?)
    input:
        score_files=rules.nested_sparse_thresholding_1kg.output,
        harmonised_geno=expand("custom_input/genotypes/{{bbid}}/chr{chr}.bed", chr=range(1,23)),
        frq = expand(rules.calculate_maf_ancestry_ext.output['afreq'], superpop='EUR', chr=range(1,23), allow_missing=True)
    output:
        profiles = "results/{bbid}/PRS_for_comparison/pt_clump/{study}/{study}.profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0],
        target_plink_chr=lambda wc, input: input['harmonised_geno'][0].replace('1.bed',''),
        ref_freq_chr=lambda wc, input: input['frq'][0][:-5],
        out_prefix=lambda wc, output: output['profiles'].replace('.profiles','')
    log:
        "logs/nested_sparse_thresholding_score_ext/{bbid}/{study}.log"
    singularity:
        config['singularity']['all']
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer/Scaled_polygenic_scorer.R "
        "--target_plink_chr {params[target_plink_chr]} "
        "--target_keep results/{wildcards[bbid]}/Ancestry_idenitfier/AllAncestry.EUR.keep "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {params[ref_freq_chr]} "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {params[out_prefix]} "
        ") &> {log} "


rule all_nested_sparse_thresholding_score_ext_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.nested_sparse_thresholding_score_ext_ref1kg.output, zip, study=studies.study_id.iloc[0], bbid=bbids)



rule sparse_thresholding_score_ext_ref1kg:
    # P+T clump sparse (non nested)
    # TODO: other ancestries
    # TODO: move input files (?)
    input:
        score_files = rules.sparse_thresholding_1kg.output,
        harmonised_geno=expand("custom_input/genotypes/{{bbid}}/chr{chr}.bed", chr=range(1,23)),
        frq = expand(rules.calculate_maf_ancestry_ext.output['afreq'], superpop='EUR', chr=range(1,23), allow_missing=True)
    output:
        profiles="results/{bbid}/PRS_for_comparison/pt_clump_nonnested/{study}/{study}.profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0],
        target_plink_chr=lambda wc, input: input['harmonised_geno'][0].replace('1.bed',''),
        ref_freq_chr=lambda wc, input: input['frq'][0][:-5],
        out_prefix=lambda wc, output: output['profiles'].replace('.profiles','')
    log:
        "logs/sparse_thresholding_score_ext/{bbid}/{study}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer/Scaled_polygenic_scorer.R "
        "--target_plink_chr {params[target_plink_chr]} "
        "--target_keep results/{wildcards[bbid]}/Ancestry_idenitfier/AllAncestry.EUR.keep "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump_nonnested/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump_nonnested/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {params[ref_freq_chr]} "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {params[out_prefix]} "
        ") &> {log} "
        

rule all_sparse_thresholding_score_ext_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.sparse_thresholding_score_ext_ref1kg.output, zip, study=studies.study_id, bbid=bbids)


rule dense_thresholding_score_ext_ref1kg:
    # TODO: outputs
    # P+T clump - dense (note: this one took long to finish according to GenoPred... might skip)
    input:
        score_files=rules.dense_thresholding_1kg.output,
        harmonised_geno=expand("custom_input/genotypes/{{bbid}}/chr{chr}.bed", chr=range(1,23)),
        frq = expand(rules.calculate_maf_ancestry_ext.output['afreq'], superpop='EUR', chr=range(1,23), allow_missing=True)
    output:
        profiles = "results/{bbid}/PRS_for_comparison/pt_clump_dense/{study}/{study}.profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0],
        target_plink_chr=lambda wc, input: input['harmonised_geno'][0].replace('1.bed',''),
        ref_freq_chr=lambda wc, input: input['frq'][0][:-5],
        out_prefix=lambda wc, output: output['profiles'].replace('.profiles','')
    log:
        "logs/dense_thresholding_score_ext/{bbid}/{study}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_dense/Scaled_polygenic_scorer_dense.R "
        "--target_plink_chr {params[target_plink_chr]} "
        "--target_keep results/{wildcards[bbid]}/Ancestry_idenitfier/AllAncestry.EUR.keep "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump_dense/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.GWAS_sumstats_clumped.txt "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump_dense/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--prsice_path {config[prsice_path]} "
        "--rscript Rscript "
        "--pheno_name {wildcards[study]} "
        "--output {params[out_prefix]} "
        ") &> {log} "

rule all_dense_thresholding_score_ext_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.dense_thresholding_score_ext_ref1kg.output, zip, study=studies.study_id, bbid=bbids)



#############
## lassosum #
#############

rule lassosum_score_ext:
    # lassosum, using biobank MAF
    # scoring
    # TODO: make this work for all target ancestries
    # TODO: move input files to different location (?)
    input:
        score_files = lambda wc: expand(rules.lassosum_prep.output, study=wc['study'], ancestry=studies.ancestry[studies.study_id == wc.study].iloc[0]),
        harmonised_geno = expand("custom_input/genotypes/{{bbid}}/chr{chr}.bed", chr=range(1,23)),
        frq = expand(rules.calculate_maf_ancestry_ext.output['afreq'], superpop='EUR', chr=range(1,23), allow_missing=True)
    output:
        profiles="results/{bbid}/PRS_for_comparison/lassosum/{study}/{study}.lassosum_profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0],
        target_plink_chr=lambda wc, input: input['harmonised_geno'][0].replace('1.bed',''),
        ref_freq_chr=lambda wc, input: input['frq'][0][:-5],
        out_prefix=lambda wc, output: output['profiles'].replace('.lassosum_profiles','')
    log:
        "logs/lassosum_score_ext/{bbid}/{study}.log"
    singularity:
        config['singularity']['all']
    threads:
        8
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_lassosum/Scaled_polygenic_scorer_lassosum.R "
        "--target_plink_chr custom_input/genotypes/{wildcards[bbid]}/chr "
        "--target_keep results/{wildcards[bbid]}/Ancestry_idenitfier/AllAncestry.EUR.keep "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/lassosum/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/lassosum/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {params[ref_freq_chr]} "
        "--pheno_name {wildcards[study]} "
        "--n_cores {threads} "
        "--plink {config[plink1_9]} "
        "--output {params[out_prefix]} "
        ") &> {log} "
        # "--target_keep custom_input/phenotypes/{wildcards[bbid]}/{params[pheno]}.txt "
        
rule all_lassosum_score_ext:
    # run rule above for all studies
    input:
        expand(rules.lassosum_score_ext.output, zip, study=studies.study_id, allow_missing=True, bbid=bbids)

    
##########
## PRScs #
##########


rule prscs_score_ext_refukbb:
    # Scoring using the pre-computed UKBB reference 
    # note: this is different from the original GenoPred
    # the scores are still normalized to 1000G (not to UKBB)!
    # however, the LD reference is the pre-computed UKBB reference
    # Note: this will crash / not run if study_ancestry is not EUR
    # TODO: make this work for all target ancestries
    # TODO: move input files to different location (?)
    input:
        score_files=rules.run_prscs_precompld_1kg_refukbb.output,
        harmonised_geno=expand("custom_input/genotypes/{{bbid}}/chr{chr}.bed", chr=range(1,23)),
        frq = expand(rules.calculate_maf_ancestry_ext.output['afreq'], superpop='EUR', chr=range(1,23), allow_missing=True)
    output:
        profiles="results/{bbid}/PRS_for_comparison/PRScs/{study}/{study}.PRScs_profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0],
        target_plink_chr=lambda wc, input: input['harmonised_geno'][0].replace('1.bed',''),
        ref_freq_chr=lambda wc, input: input['frq'][0][:-5],
        out_prefix=lambda wc, output: output['profiles'].replace('.PRScs_profiles','')
    log:
        "logs/prscs_score_ext/{bbid}/{study}.log"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_PRScs/Scaled_polygenic_scorer_PRScs.R "
        "--target_plink_chr {params[target_plink_chr]} "
        "--target_keep results/{wildcards[bbid]}/Ancestry_idenitfier/AllAncestry.EUR.keep "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/PRScs_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/PRScs_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {params[ref_freq_chr]} "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {params[out_prefix]} "
        ") &> {log} "
        

rule all_prscs_score_ext_refukbb:
    # run rule above for all studies
    input:
        expand(rules.prscs_score_ext_refukbb.output, zip, study=studies.study_id, bbid=bbids)
            

##########
## SBLUP #
##########

rule sblup_score_ext_ref1kg:
    # SBLUP
    # Scoring using the 1000G reference
    # TODO: make this work for all target ancestries
    # TODO: move input files to different location (?)
    input:
        score_files = lambda wc: expand(rules.sblup_prep.output, study=wc['study'], ancestry=studies.ancestry[studies.study_id == wc.study].iloc[0]),
        harmonised_geno=expand("custom_input/genotypes/{{bbid}}/chr{chr}.bed", chr=range(1,23)),
        frq = expand(rules.calculate_maf_ancestry_ext.output['afreq'], superpop='EUR', chr=range(1,23), allow_missing=True)
    output:
        profiles="results/{bbid}/PRS_for_comparison/SBLUP/{study}/{study}.SBLUP_profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0],
        target_plink_chr=lambda wc, input: input['harmonised_geno'][0].replace('1.bed',''),
        ref_freq_chr=lambda wc, input: input['frq'][0][:-5],
        out_prefix=lambda wc, output: output['profiles'].replace('.SBLUP_profiles','')
    log:
        "logs/sblup_score_ext_ref1kg/{bbid}/{study}.log"
    singularity:
        config['singularity']['all']
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_SBLUP/Scaled_polygenic_scorer_SBLUP.R "
        "--target_plink_chr {params[target_plink_chr]} "
        "--target_keep results/{wildcards[bbid]}/Ancestry_idenitfier/AllAncestry.EUR.keep "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/sblup/{wildcards[study]}/GWAS_sumstats_SBLUP.sblup.cojo "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/sblup/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {params[ref_freq_chr]} "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {params[out_prefix]} "
        ") &> {log} "

rule all_sblup_score_ext_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.sblup_score_ext_ref1kg.output, zip, study=studies.study_id, bbid=bbids)


###########
# SBayesR #
###########

rule sbayesr_score_ext_refukbb_robust:
    # SBayesR
    # Scoring using the pre-computed UKBB reference
    # note: this is different from the original GenoPred
    # the scores are still normalized to 1000G (not to UKBB)!
    # however, the LD reference is the pre-computed UKBB reference
    # Note: this will crash if ancestry is not EUR! 
    # TODO: make this work for all target ancestries
    # TODO: move input files to different location (?)
    input:
        rules.run_sbayesr_precompld_1kg_refukbb_robust.output,
        harmonised_geno=expand("custom_input/genotypes/{{bbid}}/chr{chr}.bed", chr=range(1,23)),
        frq=expand(rules.calculate_maf_ancestry_ext.output['afreq'], superpop='EUR', chr=range(1,23), allow_missing=True)
    output:
        profiles="results/{bbid}/PRS_for_comparison/SBayesR/{study}/{study}.sbayesr_profiles",
        log="results/{bbid}/PRS_for_comparison/SBayesR/{study}/{study}.log"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0],
        target_plink_chr=lambda wc, input: input['harmonised_geno'][0].replace('1.bed',''),
        ref_freq_chr=lambda wc, input: input['frq'][0][:-5],
        out_prefix=lambda wc, output: output['profiles'].replace('.sbayesr_profiles','')
    log:
        "logs/sbayesr_score_ext_refukbb/{bbid}/{study}.log"
    singularity:
        config['singularity']['all']
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_SBayesR/Scaled_polygenic_scorer_SBayesR.R "
        "--target_plink_chr {params[target_plink_chr]} "
        "--target_keep results/{wildcards[bbid]}/Ancestry_idenitfier/AllAncestry.EUR.keep "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{wildcards[study]}/GWAS_sumstats_SBayesR.GW.snpRes "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {params[ref_freq_chr]} "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {params[out_prefix]} "
        ") &> {log} "
        

rule all_sbayesr_score_ext_refukbb_robust:
    # run rule above for all studies
    input:
        expand(rules.sbayesr_score_ext_refukbb_robust.output, zip, study=studies.study_id, bbid=bbids)


##########
# LDpred #
##########

rule ldpred_score_ext_ref1kg:
    # LDpred
    # scoring using the 1000G reference
    # TODO: make this work for all target ancestries
    # TODO: move input files to different location (?)
    input:
        score_files = lambda wc: expand(rules.ldpred_prep.output, study=wc['study'], ancestry=studies.ancestry[studies.study_id == wc.study].iloc[0]),
        harmonised_geno=expand("custom_input/genotypes/{{bbid}}/chr{chr}.bed", chr=range(1,23)),
        frq = expand(rules.calculate_maf_ancestry_ext.output['afreq'], superpop='EUR', chr=range(1,23), allow_missing=True)
    output:
        profiles="results/{bbid}/PRS_for_comparison/LDPred/{study}/{study}.LDPred_profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0],
        target_plink_chr=lambda wc, input: input['harmonised_geno'][0].replace('1.bed',''),
        ref_freq_chr=lambda wc, input: input['frq'][0][:-5],
        out_prefix=lambda wc, output: output['profiles'].replace('.LDPred_profiles','')
    log:
        "logs/ldpred_score_ext_ref1kg/{bbid}/{study}.log"
    singularity:
        config['singularity']['all']
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_LDPred/Scaled_polygenic_scorer_LDPred.R "
        "--target_plink_chr {params[target_plink_chr]} "
        "--target_keep results/{wildcards[bbid]}/Ancestry_idenitfier/AllAncestry.EUR.keep "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/ldpred/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/ldpred/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {params[ref_freq_chr]} "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {params[out_prefix]} "
        ") &> {log} "

rule all_ldpred_score_ext_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.ldpred_score_ext_ref1kg.output, zip, study=studies.study_id, bbid=bbids)


############
# LDpred 2 #
############

#rule ldpred2_score_ext_refukbb:
#    # LDpred 2
#    # scoring using the pre-computed UKBB reference
#    # note: won't work for non-EUR
#    # note: this is different from the original GenoPred
#    # the scores are still normalized to 1000G
#    # however, the LD reference is the pre-computed UKBB reference
#    # TODO: make this work for all target ancestries
#    # TODO: move input files to different location (?)
#    input:
#        rules.run_ldpred2_precompld_1kg.output,
#        harmonised_geno=expand("custom_input/genotypes/{{bbid}}/chr{chr}.bed", chr=range(1,23)),
#        frq=expand(rules.calculate_maf_ancestry_ext.output['afreq'], superpop='EUR', chr=range(1,23), allow_missing=True)
#    output:
#        profiles = "results/{bbid}/PRS_for_comparison/LDPred2/{study}/{study}.LDPred_profiles"
#    params:
#        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
#        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0],
#        target_plink_chr=lambda wc, input: input['harmonised_geno'][0].replace('1.bed',''),
#        ref_freq_chr=lambda wc, input: input['frq'][0][:-5],
#        out_prefix=lambda wc, output: output['profiles'].replace('.LDPred_profiles','')
#    log:
#        "logs/ldpred2_score_ext_refukbb/{bbid}/{study}.log"
#    singularity:
#        config['singularity']['all']
#    threads:
#        6
#    shell:
#        "( "
#        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_LDPred2/Scaled_polygenic_scorer_LDPred2.R "
#        "--target_plink_chr {params[target_plink_chr]} "
#        "--target_keep results/{wildcards[bbid]}/Ancestry_idenitfier/AllAncestry.EUR.keep "
#        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/LDPred2_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
#        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/LDPred2_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
#        "--ref_freq_chr {params[ref_freq_chr]} "
#        "--plink {config[plink1_9]} "
#        "--n_cores {threads} "
#        "--pheno_name {wildcards[study]} "
#        "--output {params[out_prefix]} "
#        ") &> {log} "
#        
#rule all_ldpred2_score_ext_refukbb:
#    # run rule above for all studies
#    input:
#        expand(rules.ldpred2_score_ext_refukbb.output, zip, study=studies.study_id, bbid=bbids)
#

##########
# DBSLMM #
##########

rule dbslmm_score_ext_ref1kg:
    # DBSLMM
    # scoring using the 1000G reference
    input:
        score_files = lambda wc: expand(rules.dbslmm_prep.output, study=wc['study'], ancestry=studies.ancestry[studies.study_id == wc.study].iloc[0]),
        harmonised_geno=expand("custom_input/genotypes/{{bbid}}/chr{chr}.bed", chr=range(1,23)),
        frq = expand(rules.calculate_maf_ancestry_ext.output['afreq'], superpop='EUR', chr=range(1,23), allow_missing=True)
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0],
        target_plink_chr=lambda wc, input: input['harmonised_geno'][0].replace('1.bed',''),
        ref_freq_chr=lambda wc, input: input['frq'][0][:-5],
        out_prefix=lambda wc, output: output['profiles'].replace('.DBSLMM_profiles','')
    output:
        profiles="results/{bbid}/PRS_for_comparison/DBSLMM/{study}/{study}.DBSLMM_profiles"
    log:
        "logs/dbslmm_score_ext_ref1kg/{bbid}/{study}.log"
    singularity:
        config['singularity']['all']
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_DBSLMM/Scaled_polygenic_scorer_DBSLMM.R "
        "--target_plink_chr {params[target_plink_chr]} "
        "--target_keep results/{wildcards[bbid]}/Ancestry_idenitfier/AllAncestry.EUR.keep "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/dbslmm/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.dbslmm.GW.txt "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/dbslmm/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {params[ref_freq_chr]} "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {params[out_prefix]} "
        ") &> {log} "
        
rule all_dbslmm_score_ext_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.dbslmm_score_ext_ref1kg.output, zip, study=studies.study_id, bbid=bbids)


########
## ALL #
########

# rule all_score_ext:
#     # rule that runs all the rules above
#     input:
#         rules.all_sparse_thresholding_score_ext_ref1kg.input,
#         rules.all_lassosum_score_ext.input,
#         rules.all_prscs_score_ext_refukbb.input,
#         rules.all_sblup_score_ext_ref1kg.input,
#         rules.all_sbayesr_score_ext_refukbb_robust.input,
#         rules.all_ldpred_score_ext_ref1kg.input,
#         rules.all_ldpred2_score_ext_refukbb.input,
#         rules.all_dbslmm_score_ext_ref1kg.input
# 
        
        
###########################
# END: Polygenic scoring  #
###########################


############################>>
# START: Score evaluation  #>>
############################>>


#########
# P + T #
#########

# TODO

#########
# Other #
#########

#rule model_eval_ext_prep:
#    # methods are evaluated *together*
#    # phenotypes are evaluated *separately*
#    input:
#        pt_clump = rules.sparse_thresholding_score_ext_ref1kg.output['profiles'],
#        lassosum = rules.lassosum_score_ext.output['profiles'],
#        prscs = rules.prscs_score_ext_refukbb.output['profiles'],
#        sblup = rules.sblup_score_ext_ref1kg.output['profiles'],
#        sbayesr = rules.sbayesr_score_ext_refukbb_robust.output['profiles'],
#        dbslmm = rules.dbslmm_score_ext_ref1kg.output['profiles'],
#        ldpred = rules.ldpred_score_ext_ref1kg.output['profiles'],
#        ldpred2 = rules.ldpred2_score_ext_refukbb.output['profiles']
#    output:
#        predictors = 'results/{bbid}/PRS_evaluation/{study}/{study}.AllMethodComp.predictor_groups'
#    run:
#        # prepare a file with predictors, grouped by method
#        with open(output[0], 'w') as outfile:
#            outfile.write('predictors group\n')
#            for k, v in input.items():
#                outfile.write('{} {}\n'.format(v, k))
                
#rule all_model_eval_ext_prep:
#    input:
#        expand(rules.model_eval_ext_prep.output, study=studies.study_id, allow_missing=True)
#
#
#rule model_eval_ext:
#    input:
#        predictors = rules.model_eval_ext_prep.output.predictors,
#        pheno_file = lambda wc: 'custom_input/phenotypes/{bbid}/' + studies.name[studies.study_id == wc.study].iloc[0] +'.txt'
#    output:
#        assoc='results/{bbid}/PRS_evaluation/{study}/{study}.AllMethodComp.assoc.txt',
#        pred_comp='results/{bbid}/PRS_evaluation/{study}/{study}.AllMethodComp.pred_comp.txt',
#        pred_eval='results/{bbid}/PRS_evaluation/{study}/{study}.AllMethodComp.pred_eval.txt'
#    params:
#        prev = lambda wc: prevalence[ studies.name[studies.study_id == wc['study']].iloc[0] ],
#        out_prefix = lambda wc, output: output['assoc'].replace('.assoc.txt','')
#    threads:
#        8
#    log:
#        'logs/model_eval_ext/{bbid}/{study}.log'
#    singularity:
#        config['singularity']['all']
#    shell:
#        "("
#        "Rscript {config[GenoPred_dir]}/Scripts/Model_builder/Model_builder_V2.R "
#        "--pheno {input[pheno_file]} "
#        "--predictors {input[predictors]} "
#        "--n_core {threads} "
#        "--compare_predictors T "
#        "--assoc T "
#        "--outcome_pop_prev {params[prev]} "
#        "--out {params[out_prefix]} "
#        "--save_group_model T "
#        ") &> {log}"


# rule all_model_eval_ext:
#     input:
#         expand(rules.model_eval_ext.output, study=studies.study_id, allow_missing=True)
# 
        
# rule get_best_models_ext:
#     # exctract the best performing models and their performance into a neat TSV
#     input:
#         pred_eval=rules.model_eval_ext.output['pred_eval']
#     output:
#         best_models_tsv='results/{bbid}/PRS_evaluation/{study}/best_models.tsv'
#     log:
#         'logs/get_best_models_ext/{bbid}/{study}.log'
#     singularity:
#         config['singularity']['all']
#     shell:
#         "("
#         "Rscript workflow/scripts/R/get_best_models_from_pred_eval.R "
#         "--pred_eval {input[pred_eval]} "
#         "--drop TRUE "
#         ") &> {log}"
        
# rule all_get_best_models_ext:
#     input:
#         expand(rules.get_best_models_ext.output, study=studies.study_id, bbid=bbids)


# rule model_eval_custom_ext:
#     # Run Model_builder_V2.R for custom group - phenotype combinations
#     input:
#         predictors = 'custom_input/{bbid}/predictor_groups/{group}.predictor_groups',
#         pheno_file = 'custom_input/{bbid}/phenotypes/{pheno}.txt'
#     output:
#         assoc='results/{bbid}/PRS_evaluation_custom/{group}/{pheno}/AllMethodComp.assoc.txt',
#         pred_comp='results/{bbid}/PRS_evaluation_custom/{group}/{pheno}/AllMethodComp.pred_comp.txt',
#         pred_eval='results/{bbid}/PRS_evaluation_custom/{group}/{pheno}/AllMethodComp.pred_eval.txt'
#     params:
#         prev = lambda wc: prevalence[wc['pheno']],
#         out_prefix = lambda wc, output: output['assoc'].replace('.assoc.txt','')
#     threads:
#         16
#     log:
#         'logs/model_eval_custom_ext/{bbid}/{group}_{pheno}.log'
#     singularity:
#         config['singularity']['all']
#     shell:
#         "("
#         "Rscript {config[GenoPred_dir]}/Scripts/Model_builder/Model_builder_V2.R "
#         "--pheno {input[pheno_file]} "
#         "--predictors {input[predictors]} "
#         "--n_core {threads} "
#         "--compare_predictors T "
#         "--assoc T "
#         "--outcome_pop_prev {params[prev]} "
#         "--out {params[out_prefix]} "
#         "--save_group_model T "
#         ") &> {log}"
# 

# rule get_best_models_custom_ext:
#     # exctract the best performing models and their performance into a neat TSV
#     input:
#         pred_eval=rules.model_eval_custom_ext.output['pred_eval']
#     output:
#         best_models_tsv='results/{bbid}/PRS_evaluation_custom/{group}/{pheno}/best_models.tsv'
#     log:
#         'logs/get_best_models_custom/{bbid}/{group}_{pheno}.log'
#     singularity:
#         config['singularity']['all']
#     shell:
#         "("
#         "Rscript workflow/scripts/R/get_best_models_from_pred_eval.R "
#         "--pred_eval {input[pred_eval]} "
#         "--drop TRUE "
#         ") &> {log} "
#     
