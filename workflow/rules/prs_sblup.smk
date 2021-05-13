# Rules for preparing score and scale files for polygenic scoring using sblup

rule install_sblup_software:
    # Install the external software required by the sblup method
    output:
        ldsc_software="{}/ldsc.py".format(config['LDSC_dir']),
        munge_sumstats_software="{}/munge_sumstats.py".format(config['LDSC_dir'])
    shell:
        "bash install_software.sh"


rule download_sblup_ld_ref:
    # Download the LD refence panel for the given ancestry (note that it currently only supports EUR and EAS) 
    output:
        ld_ref="{}/{{ancestry}}/1.l2.ldscore.gz".format(config['LD_ref_dir'])
    shell:
        "bash download_resources.sh"


rule sblup_prep:
    # Note that LDSC requires python 2 so Snakemake will setup this environment using the given .yml file
    input: 
        ldsc_software=rules.install_sblup_software.output.ldsc_software,
        munge_sumstats_software=rules.install_sblup_software.output.munge_sumstats_software,
        ld_ref=rules.download_sblup_ld_ref.output.ld_ref,
        sumstats="{}/{{study}}.{{ancestry}}.cleaned.gz".format(config['Base_sumstats_dir'])
    output:
        "{}/sblup/{{study}}/1KGPhase3.w_hm3.{{ancestry}}.{{study}}.{{ancestry}}.scale.{{ancestry}}.scale".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/prs_sblup_{study}.{ancestry}.log"
    conda:
        "../../{}/environment.yml".format(config['LDSC_dir'])
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_SBLUP/polygenic_score_file_creator_SBLUP.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input.sumstats} "
        "--plink {config[plink1_9]} "
        "--gcta {config[gcta]} "
        "--munge_sumstats {input.munge_sumstats_software} "
        "--ldsc {input.ldsc_software} "
        "--ldsc_ref {config[LD_ref_dir]}/{wildcards[ancestry]} "
        "--hm3_snplist {config[HapMap3_snplist_dir]}/w_hm3.snplist "
        "--memory 50000 "
        "--n_cores 6 "
        "--output {output} "
        "--ref_pop_scale {config[Geno_1KG_dir]}/super_pop_keep.list "
        ") &> {log}"


rule all_sblup_prep:
    input: 
        expand("{}/sblup/{{study.study_id}}/1KGPhase3.w_hm3.{{study.ancestry}}.{{study.study_id}}.{{study.ancestry}}.scale.{{study.ancestry}}.scale".format(config['Base_sumstats_dir']), study=studies.itertuples())
