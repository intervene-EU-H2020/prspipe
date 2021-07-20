

# https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html

############################
# Harmonization to HapMap3 #
############################

rule harmonize_ukbb:
    # 1.1 UK Biobank
    input:
        rules.extract_hm3.output
    output:
        log='{}/Genotype/Harmonised/UKBB.harmonisation.chr{{chr}}.log'.format(config['UKBB_output']),
        bim='{}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr{{chr}}.bim'.format(config['UKBB_output']),
        bed='{}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr{{chr}}.bed'.format(config['UKBB_output']),
        fam='{}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr{{chr}}.fam'.format(config['UKBB_output'])
    log:
        "logs/harmonize_ukbb_chr{chr}.log"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Harmonisation_of_UKBB/Harmonisation_of_UKBB.R "
        "--chr {wildcards[chr]} "
        "--input_dir {config[UKBB_original]} "
        "--target_fam {config[UKBB_fam]} "
        "--output_dir {config[UKBB_output]}/Genotype/Harmonised "
        "--reference_dir {config[Geno_1KG_dir]} "
        "--qctool2 {config[qctool2]} "
        "--plink {config[plink1_9]} "
        ") &> {log} "


####################
# Ancestry scoring #
####################

rule ancestry_scoring_ukbb:
    # 2.1 UK Biobank
    # recommended 15G of RAM
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23))
    output:
        # touch("{}/Projected_PCs/Ancestry_idenitfier/test.ok".format(config['UKBB_output'])),
        eigenvec_sp=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.{superpop}.eigenvec', ukbb_output=config['UKBB_output'], superpop=config['1kg_superpop']),
        keep_sp=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.{superpop}.keep', ukbb_output=config['UKBB_output'], superpop=config['1kg_superpop']),
        scale_sp=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.{superpop}.scale', ukbb_output=config['UKBB_output'], superpop=config['1kg_superpop']),
        eigenvec=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.eigenvec', ukbb_output=config['UKBB_output']),
        eigenvec_var=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.eigenvec.var', ukbb_output=config['UKBB_output']),
        log=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.log', ukbb_output=config['UKBB_output']),
        model_pred_keep=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.model_pred.{superpop}.keep', ukbb_output=config['UKBB_output'], superpop=config['1kg_superpop']),
        model_pred=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.model_pred', ukbb_output=config['UKBB_output']),
        pop_model=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.pop_model.rds', ukbb_output=config['UKBB_output']),
        model_details=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.pop_model_prediction_details.txt', ukbb_output=config['UKBB_output']),
        scale=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.scale', ukbb_output=config['UKBB_output'])
    log:
        "logs/ancestry_scoring_ukbb.log"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Ancestry_identifier/Ancestry_identifier.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--n_pcs 100 "
        "--plink {config[plink1_9]} "
        "--plink2 {config[plink2]} "
        "--output {config[UKBB_output]}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--pop_data {config[Geno_1KG_dir]}/integrated_call_samples_v3.20130502.ALL.panel_small "
        ") &> {log} "
        
        
rule ancestry_scoring_ukbb_stringent:
    # get ancestry keep files using a more stringent cutoff of 0.995
    input:
        rules.ancestry_scoring_ukbb.output['model_pred']
    output:
        keep=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.{superpop}.keep', ukbb_output=config['UKBB_output'], superpop=config['1kg_superpop'])
    params:
        ancestries=" ".join(sorted(config['1kg_superpop']))
    log:
        'logs/ancestry_scoring_ukbb_stringent.log'
    shell:
        "("
        "bash workflow/scripts/bash/extract_ancestry.sh {input} \"{params[ancestries]}\" {config[UKBB_output]}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry "
        ") &> {log}"
        
        
rule calculate_maf_ancestry_ukbb:
    # calculate allele frequencies for different superpopulations
    # output of this rule could be used if we wanted to use freq-files for other ancestries than EUR
    input:
        rules.ancestry_scoring_ukbb_stringent.output,
        rules.harmonize_ukbb.output
    output:
        afreq_gz='{}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.{{superpop}}.chr{{chr}}.frq.gz'.format(config['UKBB_output'])
    log:
        'logs/calculate_maf_ancestry_ukbb_{superpop}_{chr}.log'
    shell:
        '('
        '{config[plink1_9]} '
        '--bfile {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr{wildcards[chr]} '
        '--freq '
        '--keep {config[UKBB_output]}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.{wildcards[superpop]}.keep '
        '--out {config[UKBB_output]}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.{wildcards[superpop]}.chr{wildcards[chr]}; '
        'gzip {config[UKBB_output]}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.{wildcards[superpop]}.chr{wildcards[chr]}.frq; '
        ') &> {log} '
    

rule all_calculate_maf_ancestry_ukbb:
    # runs rule above for all ancestries and chromosomes
    input:
        expand(rules.calculate_maf_ancestry_ukbb.output, chr=range(1,23), superpop=config['1kg_superpop'])
    output:
        touch('{}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.all_afreq.ok'.format(config['UKBB_output']))
        

######################
# UKBB 10K reference #
######################

rule export_ukbb10k_keep:
    # generate keep-file for EUR UKBB 10K reference
    # https://opain.github.io/GenoPred/Pipeline_prep_withUKBB_ref.html
    # 1.1
    input:
        expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        rules.ancestry_scoring_ukbb.output
    output:
        config['UKBB_output'] + '/UKBB_ref/keep_files/UKBB_noPheno_EUR_10K.keep'
    log:
        'logs/export_ukbb10k_keep/EUR.log'
    shell:
        "("
        "Rscript workflow/scripts/R/ukbb/export_ukbb10k_keep.R"
        ") &> {log}"


rule extract_ukbb10k:
    # extract genotypes of EUR UKBB 10K subset and calculate MAF
    # https://opain.github.io/GenoPred/Pipeline_prep_withUKBB_ref.html
    # 2 
    input:
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        keep_files=rules.export_ukbb10k_keep.output
    output:
        bed=expand(config['UKBB_output'] + '/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr{chr}.bed', chr=range(1,23)),
        bim=expand(config['UKBB_output'] + '/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr{chr}.bim', chr=range(1,23)),
        fam=expand(config['UKBB_output'] + '/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr{chr}.fam', chr=range(1,23)),
        frq=expand(config['UKBB_output'] + '/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr{chr}.frq', chr=range(1,23)),
        bed_gw = config['UKBB_output'] + '/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.GW.bed',
        bim_gw = config['UKBB_output'] + '/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.GW.bim',
        fam_gw = config['UKBB_output'] + '/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.GW.fam'
    log:
        "logs/extract_ukbb10k/EUR.log"
    shell:
        "("
        "for chr in $(seq 1 22); do "
        "{config[plink1_9]} "
        "--bfile {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr${{chr}} "
        "--make-bed "
        "--geno 0.02 "
        "--keep {config[UKBB_output]}/UKBB_ref/keep_files/UKBB_noPheno_EUR_10K.keep "
        "--out {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr${{chr}}; "
        " done; "
        "for chr in $(seq 1 22); do "
        "{config[plink1_9]} "
        "--bfile {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr${{chr}} "
        "--freq "
        "--out {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr${{chr}}; "
        "done; "
        "ls {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr*.bed | sed -e 's/.bed//g' > {config[UKBB_output]}/UKBB_ref/genotype/merge_list.txt; "
        "{config[plink1_9]} "
        "--merge-list {config[UKBB_output]}/UKBB_ref/genotype/merge_list.txt "
        "--make-bed "
        "--out {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.GW; "
        "rm {config[UKBB_output]}/UKBB_ref/genotype/merge_list.txt; "
        ") &> {log} "

        
#################################>>
# START: Polygenic scoring UKBB #>>
#################################>>

rule ukbb_prs_input:
    # rule that checks if the phenotype input files are present
    input:
        phenotypes = expand('{}/Phenotype/PRS_comp_subset/UKBB.{{pheno}}.txt'.format(config['UKBB_output']), pheno=studies.name)


##########################
# Pruning & Thresholding #
##########################

       
rule nested_sparse_thresholding_score_ukbb_ref1kg:
    # P+T clump - nested sparse 
    # TODO: outputs
    input:
        score_files=rules.nested_sparse_thresholding_1kg.output,
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        frq = rules.extract_ukbb10k.output.frq
    output:
        profiles = config['UKBB_output'] + "/PRS_for_comparison/1KG_ref/pt_clump/{study}/UKBB.subset.w_hm3.{study}.profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0]
    log:
        "logs/nested_sparse_thresholding_score_ukbb/{study}.log"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer/Scaled_polygenic_scorer.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--target_keep {config[UKBB_output]}/Phenotype/PRS_comp_subset/UKBB.{params[pheno]}.subsample.txt "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {config[UKBB_output]}/PRS_for_comparison/1KG_ref/pt_clump/{wildcards[study]}/UKBB.subset.w_hm3.{wildcards[study]} "
        ") &> {log} "
        # force UKBB MAF not, 1KG:
        # "--ref_freq_chr {config[Geno_1KG_dir]}/freq_files/{params[study_ancestry]}/1KGPhase3.w_hm3.{params[study_ancestry]}.chr "


rule all_nested_sparse_thresholding_score_ukbb_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.nested_sparse_thresholding_score_ukbb_ref1kg.output, zip, study=studies.study_id.iloc[0] )

        
rule dense_thresholding_score_ukbb_ref1kg:
    # TODO: outputs
    # P+T clump - dense (note: this one took long to finish according to GenoPred... might skip)
    input:
        score_files=rules.dense_thresholding_1kg.output,
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23))
    output:
        profiles = config['UKBB_output'] + "/PRS_for_comparison/1KG_ref/pt_clump_dense/{study}/UKBB.subset.w_hm3.{study}.profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0]
    log:
        "logs/dense_thresholding_score_ukbb/{study}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_dense/Scaled_polygenic_scorer_dense.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--target_keep {config[UKBB_output]}/Phenotype/PRS_comp_subset/UKBB.{params[pheno]}.subsample.txt "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump_dense/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.GWAS_sumstats_clumped.txt "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump_dense/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--prsice_path {config[prsice_path]} "
        "--rscript Rscript "
        "--pheno_name {wildcards[study]} "
        "--output {config[UKBB_output]}/PRS_for_comparison/1KG_ref/pt_clump_dense/{wildcards[study]}/UKBB.subset.w_hm3.{wildcards[study]} "
        ") &> {log} "

rule all_dense_thresholding_score_ukbb_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.dense_thresholding_score_ukbb_ref1kg.output, zip, study=studies.study_id.iloc[0] )


rule sparse_thresholding_score_ukbb_ref1kg:
    # P+T clump sparse (non nested)
    # TODO: outputs
    input:
        score_files = rules.sparse_thresholding_1kg.output,
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        frq = rules.extract_ukbb10k.output.frq
    output:
        profiles=config['UKBB_output'] + "/PRS_for_comparison/1KG_ref/pt_clump_nonnested/{study}/UKBB.subset.w_hm3.{study}.profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0]
    log:
        "logs/sparse_thresholding_score_ukbb/{study}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer/Scaled_polygenic_scorer.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--target_keep {config[UKBB_output]}/Phenotype/PRS_comp_subset/UKBB.{params[pheno]}.subsample.txt "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump_nonnested/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/pt_clump_nonnested/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {config[UKBB_output]}/PRS_for_comparison/1KG_ref/pt_clump_nonnested/{wildcards[study]}/UKBB.subset.w_hm3.{wildcards[study]} "
        ") &> {log} "
        

rule all_sparse_thresholding_score_ukbb_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.sparse_thresholding_score_ukbb_ref1kg.output, zip, study=studies.study_id.iloc[0] )


############
# lassosum #
############

rule lassosum_score_ukbb:
    # lassosum, using 1000G reference MAF
    # scoring
    # TODO: outputs
    input:
        score_files = lambda wc: expand(rules.lassosum_prep.output, study=wc['study'], ancestry=studies.ancestry[studies.study_id == wc.study].iloc[0]),
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        frq = rules.extract_ukbb10k.output.frq
    output:
        profiles=config['UKBB_output'] + "/PRS_for_comparison/1KG_ref/lassosum/{study}/UKBB.subset.w_hm3.{study}.lassosum_profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0]
    log:
        "logs/lassosum_score_ukbb/{study}.log"
    threads:
        10
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_lassosum/Scaled_polygenic_scorer_lassosum.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--target_keep {config[UKBB_output]}/Phenotype/PRS_comp_subset/UKBB.{params[pheno]}.subsample.txt "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/lassosum/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/lassosum/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr "
        "--pheno_name {wildcards[study]} "
        "--n_cores {threads} "
        "--plink {config[plink1_9]} "
        "--output {config[UKBB_output]}/PRS_for_comparison/1KG_ref/lassosum/{wildcards[study]}/UKBB.subset.w_hm3.{wildcards[study]} "
        ") &> {log} "
        
rule all_lassosum_score_ukbb:
    # run rule above for all studies
    input:
        expand(rules.lassosum_score_ukbb.output, zip, study=studies.study_id.iloc[0])

        
        
#########
# PRScs #
#########

rule prscs_score_ukbb_ref1kg:
    # PRScs 
    # Scoring using the pre-computed 1000G reference
    input:
        score_files=rules.run_prscs_precompld_1kg.output,
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        frq = rules.extract_ukbb10k.output.frq
    output:
        profiles=config['UKBB_output'] + "/PRS_for_comparison/1KG_ref/PRScs/{study}/UKBB.subset.w_hm3.{study}.PRScs_profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0]
    log:
        "logs/prscs_score_ukbb/{study}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_PRScs/Scaled_polygenic_scorer_PRScs.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--target_keep {config[UKBB_output]}/Phenotype/PRS_comp_subset/UKBB.{params[pheno]}.subsample.txt "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/PRScs_precompld/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/PRScs_precompld/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {config[UKBB_output]}/PRS_for_comparison/1KG_ref/PRScs/{wildcards[study]}/UKBB.subset.w_hm3.{wildcards[study]} "
        ") &> {log} "

rule all_prscs_score_ukbb_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.prscs_score_ukbb_ref1kg.output, zip, study=studies.study_id.iloc[0] )

rule prscs_score_ukbb_refukbb:
    # Scoring using the pre-computed UKBB reference 
    # note: this is different from the original GenoPred
    # the scores are still normalized to 1000G (not to UKBB)!
    # however, the LD reference is the pre-computed UKBB reference + MAF calculated on the 10K subset
    # Note: this will crash / not run if study_ancestry is not EUR
    input:
        score_files=rules.run_prscs_precompld_1kg_refukbb.output,
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        frq = rules.extract_ukbb10k.output.frq
    output:
        profiles=config['UKBB_output'] + "/PRS_for_comparison/UKBB_ref/PRScs/{study}/UKBB.subset.w_hm3.{study}.PRScs_profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0]
    log:
        "logs/prscs_score_ukbb/{study}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_PRScs/Scaled_polygenic_scorer_PRScs.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--target_keep {config[UKBB_output]}/Phenotype/PRS_comp_subset/UKBB.{params[pheno]}.subsample.txt "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/PRScs_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/PRScs_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {config[UKBB_output]}/PRS_for_comparison/UKBB_ref/PRScs/{wildcards[study]}/UKBB.subset.w_hm3.{wildcards[study]} "
        ") &> {log} "
        

rule all_prscs_score_ukbb_refukbb:
    # run rule above for all studies
    input:
        expand(rules.prscs_score_ukbb_refukbb.output, zip, study=studies.study_id.iloc[0] )

#########
# SBLUP #
#########

rule sblup_score_ukbb_ref1kg:
    # SBLUP
    # Scoring using the 1000G reference
    input:
        score_files = lambda wc: expand(rules.sblup_prep.output, study=wc['study'], ancestry=studies.ancestry[studies.study_id == wc.study].iloc[0]),
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        frq = rules.extract_ukbb10k.output.frq
    output:
        profiles=config['UKBB_output'] + "/PRS_for_comparison/1KG_ref/SBLUP/{study}/UKBB.subset.w_hm3.{study}.SBLUP_profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0]
    log:
        "logs/sblup_score_ukbb_ref1kg/{study}.log"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_SBLUP/Scaled_polygenic_scorer_SBLUP.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--target_keep {config[UKBB_output]}/Phenotype/PRS_comp_subset/UKBB.{params[pheno]}.subsample.txt "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/sblup/{wildcards[study]}/GWAS_sumstats_SBLUP.sblup.cojo "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/sblup/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {config[UKBB_output]}/PRS_for_comparison/1KG_ref/SBLUP/{wildcards[study]}/UKBB.subset.w_hm3.{wildcards[study]} "
        ") &> {log} "
        

rule all_sblup_score_ukbb_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.sblup_score_ukbb_ref1kg.output, zip, study=studies.study_id.iloc[0] )


###########
# SBayesR #
###########

# The original GenoPred contains all types of different settings for SBayesR
# we won't explore all of these there, but rather keep the default.

rule sbayesr_score_ukbb_refukbb_robust:
    # SBayesR
    # Scoring using the pre-computed UKBB reference
    # note: this is different from the original GenoPred
    # the scores are still normalized to 1000G (not to UKBB)!
    # however, the LD reference is the pre-computed UKBB reference + MAF calculated on the 10K subset
    # Note: this will crash if ancestry is not EUR! -> TODO: implement for other ancestries
    input:
        rules.run_sbayesr_precompld_1kg_refukbb_robust.output,
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        frq=rules.extract_ukbb10k.output.frq
    output:
        profiles=config['UKBB_output'] + "/PRS_for_comparison/UKBB_ref/SBayesR/{study}/UKBB.subset.w_hm3.{study}.sbayesr_profiles",
        log=config['UKBB_output'] + "/PRS_for_comparison/UKBB_ref/SBayesR/{study}/UKBB.subset.w_hm3.{study}.log"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0]
    log:
        "logs/sbayesr_score_ukbb_refukbb/{study}.log"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_SBayesR/Scaled_polygenic_scorer_SBayesR.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--target_keep {config[UKBB_output]}/Phenotype/PRS_comp_subset/UKBB.{params[pheno]}.subsample.txt "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{wildcards[study]}/GWAS_sumstats_SBayesR.GW.snpRes "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {config[UKBB_output]}/PRS_for_comparison/UKBB_ref/SBayesR/{wildcards[study]}/UKBB.subset.w_hm3.{wildcards[study]} "
        ") &> {log} "
        # if we wanted to use 1000G reference minor allele frequencies:
        # "--ref_freq_chr {config[Geno_1KG_dir]}/freq_files/{params[study_ancestry]}/1KGPhase3.w_hm3.{params[study_ancestry]}.chr "
        

rule all_sbayesr_score_ukbb_refukbb_robust:
    # run rule above for all studies
    input:
        expand(rules.sbayesr_score_ukbb_refukbb_robust.output, zip, study=studies.study_id.iloc[0] )


##########
# LDpred #
##########


rule ldpred_score_ukbb_ref1kg:
    # LDpred
    # scoring using the 1000G reference
    input:
        score_files = lambda wc: expand(rules.ldpred_prep.output, study=wc['study'], ancestry=studies.ancestry[studies.study_id == wc.study].iloc[0]),
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        frq = rules.extract_ukbb10k.output.frq
    output:
        profiles=config['UKBB_output'] + "/PRS_for_comparison/1KG_ref/LDPred/{study}/UKBB.subset.w_hm3.{study}.LDPred_profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0]
    log:
        "logs/ldpred_score_ukbb_ref1kg/{study}.log"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_LDPred/Scaled_polygenic_scorer_LDPred.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--target_keep {config[UKBB_output]}/Phenotype/PRS_comp_subset/UKBB.{params[pheno]}.subsample.txt "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/ldpred/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/ldpred/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {config[UKBB_output]}/PRS_for_comparison/1KG_ref/LDPred/{wildcards[study]}/UKBB.subset.w_hm3.{wildcards[study]} "
        ") &> {log} "

rule all_ldpred_score_ukbb_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.ldpred_score_ukbb_ref1kg.output, zip, study=studies.study_id.iloc[0] )


############
# LDpred 2 #
############

rule ldpred2_score_ukbb_refukbb:
    # LDpred 2
    # scoring using the pre-computed UKBB reference
    # note: won't work for non-EUR
    # note: this is different from the original GenoPred
    # the scores are still normalized to 1000G
    # however, the LD reference is the pre-computed UKBB reference + MAF calculated on the 10K subset
    input:
        rules.run_ldpred2_precompld_1kg.output,
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        frq=rules.extract_ukbb10k.output.frq
    output:
        profiles = config['UKBB_output'] + "/PRS_for_comparison/UKBB_ref/LDPred2/{study}/UKBB.subset.w_hm3.{study}.LDPred_profiles"
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0]
    log:
        "logs/ldpred2_score_ukbb_refukbb/{study}.log"
    threads:
        6
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_LDPred2/Scaled_polygenic_scorer_LDPred2.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--target_keep {config[UKBB_output]}/Phenotype/PRS_comp_subset/UKBB.{params[pheno]}.subsample.txt "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/LDPred2_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/LDPred2_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr "
        "--plink {config[plink1_9]} "
        "--n_cores {threads} "
        "--pheno_name {wildcards[study]} "
        "--output {config[UKBB_output]}/PRS_for_comparison/UKBB_ref/LDPred2/{wildcards[study]}/UKBB.subset.w_hm3.{wildcards[study]} "
        ") &> {log} "
        
rule all_ldpred2_score_ukbb_refukbb:
    # run rule above for all studies
    input:
        expand(rules.ldpred2_score_ukbb_refukbb.output, zip, study=studies.study_id.iloc[0] )


##########
# DBSLMM #
##########

rule dbslmm_score_ukbb_ref1kg:
    # DBSLMM
    # scoring using the 1000G reference
    input:
        score_files = lambda wc: expand(rules.dbslmm_prep.output, study=wc['study'], ancestry=studies.ancestry[studies.study_id == wc.study].iloc[0]),
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23)),
        frq = rules.extract_ukbb10k.output.frq
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        pheno=lambda wc: studies.name[studies.study_id == wc.study ].iloc[0]
    output:
        profiles=config['UKBB_output'] + "/PRS_for_comparison/1KG_ref/DBSLMM/{study}/UKBB.subset.w_hm3.{study}.DBSLMM_profiles"
    log:
        "logs/dbslmm_score_ukbb_ref1kg/{study}.log"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer_DBSLMM/Scaled_polygenic_scorer_DBSLMM.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--target_keep {config[UKBB_output]}/Phenotype/PRS_comp_subset/UKBB.{params[pheno]}.subsample.txt "
        "--ref_score {config[Geno_1KG_dir]}/Score_files_for_polygenic/dbslmm/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.dbslmm.GW.txt "
        "--ref_scale {config[Geno_1KG_dir]}/Score_files_for_polygenic/dbslmm/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.{params[study_ancestry]}.scale "
        "--ref_freq_chr {config[UKBB_output]}/UKBB_ref/genotype/UKBB.noPheno.EUR.10K.chr "
        "--plink {config[plink1_9]} "
        "--pheno_name {wildcards[study]} "
        "--output {config[UKBB_output]}/PRS_for_comparison/1KG_ref/DBSLMM/{wildcards[study]}/UKBB.subset.w_hm3.{wildcards[study]} "
        ") &> {log} "
        
rule all_dbslmm_score_ukbb_ref1kg:
    # run rule above for all studies
    input:
        expand(rules.dbslmm_score_ukbb_ref1kg.output, zip, study=studies.study_id.iloc[0] )

#######
# ALL #
#######

rule all_score_ukbb:
    # rule that runs all the rules above
    input:
        rules.all_nested_sparse_thresholding_score_ukbb_ref1kg.input,
        rules.all_dense_thresholding_score_ukbb_ref1kg.input,
        rules.all_sparse_thresholding_score_ukbb_ref1kg.input,
        rules.all_lassosum_score_ukbb.input,
        rules.all_prscs_score_ukbb_ref1kg.input,
        rules.all_prscs_score_ukbb_refukbb.input,
        rules.all_sblup_score_ukbb_ref1kg.input,
        rules.all_sbayesr_score_ukbb_refukbb_robust.input,
        rules.all_ldpred_score_ukbb_ref1kg.input,
        rules.all_ldpred2_score_ukbb_refukbb.input,
        rules.all_dbslmm_score_ukbb_ref1kg.input
        
        
###############################
# END: Polygenic scoring UKBB #
###############################


################################>>
# START: Score evaluation UKBB #>>
################################>>


#########
# P + T #
#########

# TODO

#########
# Other #
#########

rule model_eval_prep:
    # methods are evaluated *together*
    # phenotypes are evaluated *separately*
    input:
        pt_clump = rules.sparse_thresholding_score_ukbb_ref1kg.output['profiles'],
        lassosum = rules.lassosum_score_ukbb.output['profiles'],
        prscs = rules.prscs_score_ukbb_refukbb.output['profiles'],
        sblup = rules.sblup_score_ukbb_ref1kg.output['profiles'],
        sbayesr = rules.sbayesr_score_ukbb_refukbb_robust.output['profiles'],
        dbslmm = rules.dbslmm_score_ukbb_ref1kg.output['profiles'],
        ldpred = rules.ldpred_score_ukbb_ref1kg.output['profiles'],
        ldpred2 = rules.ldpred2_score_ukbb_refukbb.output['profiles']
    output:
        predictors = config['UKBB_output'] + '/PRS_for_comparison/evaluation/{study}/Association_withPRS/UKBB.w_hm3.{study}.AllMethodComp.predictor_groups'
    run:
        # prepare a file with predictors, grouped by method
        with open(output[0], 'w') as outfile:
            outfile.write('predictors group\n')
            for k, v in input.items():
                outfile.write('{} {}\n'.format(v, k))
                
rule all_model_eval_prep:
    input:
        expand(rules.model_eval_prep.output, study=studies.study_id.iloc[0])


rule model_eval:
    input:
        predictors = rules.model_eval_prep.output.predictors,
        pheno_file = lambda wc: config['UKBB_output'] + '/Phenotype/PRS_comp_subset/UKBB.' + studies.name[studies.study_id == wc.study].iloc[0] +'.subsample.txt'
    output:
        # TODO: adjust outputs
        touch(config['UKBB_output'] + '/PRS_for_comparison/evaluation/{study}/Association_withPRS/UKBB.w_hm3.{study}.AllMethodComp.ok')
    params:
        prev = lambda wc: prevalence[ studies.name[studies.study_id == wc['study']].iloc[0] ]
    threads:
        8
    log:
        'logs/model_eval/{study}.log'
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/Model_builder/Model_builder_V2.R "
        "--pheno {input[pheno_file]} "
        "--predictors {input[predictors]} "
        "--n_core {threads} "
        "--compare_predictors T "
        "--assoc T "
        "--outcome_pop_prev {params[prev]} "
        "--out {config[UKBB_output]}/PRS_for_comparison/evaluation/{wildcards[study]}/Association_withPRS/UKBB.w_hm3.{wildcards[study]}.AllMethodComp "
        ") &> {log}"

rule all_model_eval:
    input:
        expand(rules.model_eval.output, study=studies.study_id)
