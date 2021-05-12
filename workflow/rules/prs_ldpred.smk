# Rules for preparing score and scale files for polygenic scoring using ldpred

# NOTE you have to run the rules using the --conda-not-block-search-path-envvars option because of the pythonpath adjustment

rule install_ldpred_software:
    # Note: There's a bug in the latest version of ldpred (1.0.11), see https://github.com/bvilhjal/ldpred/issues/151
    # To avoid this issue you need to use v1.0.10 (specifically, commit 77084f1196239ab42c92492af85128c1c3d0d0c1)
    output:
        ldpred_software="{}/ldpred/run.py".format(config['LDPRED_dir'])
    shell:
        "bash install_software.sh"


rule ldpred_prep:
    # Note 1: the ldpred software requires Python dependencies, specified in the ldpred.yaml conda environment
    # Note 2: because the R script uses the ldpred software scripts (instead of the Python package), the script exports the module to the pythonpath
    input:
        ldpred_software=rules.install_ldpred_software.output.ldpred_software,
        sumstats="{}/{{study}}.{{ancestry}}.cleaned.gz".format(config['Base_sumstats_dir'])
    output:
        "{}/ldpred/{{study}}/1KGPhase3.w_hm3.{{ancestry}}.{{study}}.{{ancestry}}.scale".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/prs_ldpred_{study}.{ancestry}.log"
    conda:
        "../envs/ldpred.yaml"
    shell:
        "export PYTHONPATH=$PYTHONPATH:$PWD/workflow/scripts/ldpred && ("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_LDPred/polygenic_score_file_creator_LDPred.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input.sumstats} "
        "--plink {config[plink1_9]} "
        "--memory 20000 "
        "--n_cores 1 "
        "--ldpred \"python -m ldpred\" "
        "--output {output} "
        "--ref_pop_scale {config[Geno_1KG_dir]}/super_pop_keep.list "
        ") &> {log}"


rule all_ldpred_prep:
    input: 
        expand("{}/ldpred/{{study.study_id}}/1KGPhase3.w_hm3.{{study.ancestry}}.{{study.study_id}}.{{study.ancestry}}.scale".format(config['Base_sumstats_dir']), study=studies.itertuples())
