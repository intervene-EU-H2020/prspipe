# Rules for preparing score and scale files for polygenic scoring using lassosum

# TODO the output file currently there is a bit of a hack to get it to run... but doesn't accurately represent all the outputs

rule lassosum_prep:
    input: 
        "{}/{{study}}.{{ancestry}}.cleaned.gz".format(config['Base_sumstats_dir'])
    output:
        "{}/lassosum/{{study}}/1KGPhase3.w_hm3.{{ancestry}}.{{study}}.{{ancestry}}.scale".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/prs_lassosum_{study}.{ancestry}.log"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_lassosum/polygenic_score_file_creator_lassosum.R "
        "--ref_plink_gw {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input} "
        "--output {config[Base_sumstats_dir]}/lassosum/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[ancestry]}.{wildcards[study]} "
        "--plink {config[plink1_9]} "
        "--ref_pop_scale {config[Geno_1KG_dir]}/super_pop_keep.list "
        ") &> {log}"


rule all_lassosum_prep:
    input: 
        expand("{}/lassosum/{{study.study_id}}/1KGPhase3.w_hm3.{{study.ancestry}}.{{study.study_id}}.{{study.ancestry}}.scale".format(config['Base_sumstats_dir']), study=studies.itertuples())
