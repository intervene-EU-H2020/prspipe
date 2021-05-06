# Rules for preparing score and scale files for polygenic scoring using ldpred

rule install_ldpred_software:
    # Check if ldpred software has been downloaded from https://github.com/bvilhjal/ldpred
    output:
        "workflow/scripts/ldpred/ldpred/run.py"
    shell:
        "bash install_software.sh"


rule ldpred_prep:
    # Note that the ldpred software requires Python dependencies, specified in the ldpred.yaml conda environment
    input:
        ldpred=rules.install_ldpred_software.output,
        sumstats="{}/{{study}}.{{ancestry}}.cleaned.gz".format(config['Base_sumstats_dir'])
    output:
        "{}/ldpred/{{study}}/1KGPhase3.w_hm3.{{ancestry}}.{{study}}".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/prs_ldpred_{study}.{ancestry}.log"
    conda:
        "../envs/ldpred.yaml"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_LDPred/polygenic_score_file_creator_LDPred.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input.sumstats} "
        "--plink {config[plink1_9]} "
        "--memory 20000 "
        "--n_cores 1 "
        "--ldpred ./{input.ldpred} "
        "--output {output} "
        "--ref_pop_scale {config['Geno_1KG_dir']}/super_pop_keep.list "
        ") &> {log}"


rule all_ldpred_prep:
    input: 
        expand("{}/ldpred/{{study.study_id}}/1KGPhase3.w_hm3.{{study.ancestry}}.{{study.study_id}}".format(config['Base_sumstats_dir']), study=studies.itertuples())
