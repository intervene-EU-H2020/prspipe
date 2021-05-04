# Rules for preparing score and scale files for polygenic scoring using lassosum

# TODO log error says it can't find input summary stats files?

rule lassosum_prep:
    input: 
        "{}/{{study}}.{{ancestry}}.cleaned.gz".format(config['Base_sumstats_dir'])
    output:
        "{}/lassosum/1KGPhase3.w_hm3.{{ancestry}}.{{study}}".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/prs_lassosum_{study}.{ancestry}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_lassosum/polygenic_score_file_creator_lassosum.R "
        "--ref_plink_gw {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input} "
        "--output {output} "
        "--plink {config[plink1_9]} "
        "--ref_pop_scale {config[Geno_1KG_dir]}/super_pop_keep.list "
        ") &> {log}"


rule all_lassosum_prep:
    input: 
        expand("{}/lassosum/1KGPhase3.w_hm3.{{study.ancestry}}.{{study.study_id}}".format(config['Base_sumstats_dir']), study=studies.itertuples())
