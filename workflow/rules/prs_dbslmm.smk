# Rules for preparing score and scale files for polygenic scoring using dbslmm

# TODO need to set ldsc, ldsc_ref, ld_blocks, r_script (?), dbslmm, sample_prev, pop_prev in config or as parameters
# TODO update install software .sh file to download dbslmm etc.

rule dbslmm_prep:
    input: 
        "{}/{{study}}.{{ancestry}}.cleaned.gz".format(config['Base_sumstats_dir'])
    output:
        "{}/dbslmm/1KGPhase3.w_hm3.{{ancestry}}.{{study}}".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/prs_dbslmm_{study}.{ancestry}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_DBSLMM/polygenic_score_file_creator_DBSLMM.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input} "
        "--plink {config[plink1_9]} "
        "--memory 5000 "
        "--ld_blocks {config[ld_blocks]} "
        "--rscript {config[r_script]} "
        "--dbslmm {config[dbslmm]} "
        "--munge_sumstats {config[GenoPred_dir]}/ldsc "
        "--ldsc {config[ldsc]} "
        "--ldsc_ref {config[ldsc_ref]} "
        "--hm3_snplist {config[HapMap3_snplist_dir]}/w_hm3.snplist "
        "--sample_prev {config[sample_prev]} "
        "--pop_prev {config[pop_prev]} "
        "--output {output} "
        "--ref_pop_scale {config['Geno_1KG_dir']}/super_pop_keep.list "
        ") &> {log}"


rule all_dbslmm_prep:
    input: 
        expand("{}/dbslmm/1KGPhase3.w_hm3.{{study.ancestry}}.{{study.study_id}}".format(config['Base_sumstats_dir']), study=studies.itertuples())
