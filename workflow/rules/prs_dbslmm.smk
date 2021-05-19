# Rules for preparing score and scale files for polygenic scoring using dbslmm

rule install_dbslmm_software:
    # Install the external software required by the dbslmm method
    # Note that we use our own fork of DBSLMM (had to make changes for compatibility with GenoPred scripts)
    output:
        ldsc_software="{}/ldsc.py".format(config['LDSC_dir']),
        munge_sumstats_software="{}/munge_sumstats.py".format(config['LDSC_dir']),
        dbslmm_software="{}/software/DBSLMM.R".format(config['DBSLMM_dir'])
    shell:
        "bash install_software.sh"


rule download_dbslmm_ld_block:
    # Download the LD block data for the given ancestry (note that it currently only supports AFR, ASN and EUR) 
    output:
        ld_block="{}/{{ancestry}}/fourier_ls-all.bed".format(config['LD_block_dir'])
    shell:
        "bash download_resources.sh"


rule dbslmm_prep:
    # Implements the dbslmm method
    # Note that LDSC requires python 2 so Snakemake will setup this environment using the given .yml file
    # Uses the precomputed LD ref (based on 1000G) by default - 
    # If you want to compute the LD ref from scratch, replace ld_ref with expand("{}/sblup_dbslmm/1000G/fromscratch/{{study.ancestry}}/{{chr_id}}.l2.ldscore.gz".format(config['LD_ref_dir']), chr_id=range(1,23), study=studies.itertuples())
    input: 
        ldsc_software=rules.install_dbslmm_software.output.ldsc_software,
        munge_sumstats_software=rules.install_dbslmm_software.output.munge_sumstats_software,
        dbslmm_software=rules.install_dbslmm_software.output.dbslmm_software,
        ld_ref=expand("{}/sblup_dbslmm/1000G/precomputed/{{study.ancestry}}/{{chr_id}}.l2.ldscore.gz".format(config['LD_ref_dir']), chr_id=range(1,23), study=studies.itertuples()),
        ld_block=rules.download_dbslmm_ld_block.output.ld_block,
        sumstats="{}/{{study}}.{{ancestry}}.cleaned.gz".format(config['Base_sumstats_dir'])
    output:
        "{}/dbslmm/{{study}}/1KGPhase3.w_hm3.{{ancestry}}.{{study}}.{{ancestry}}.scale".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/prs_dbslmm_{study}.{ancestry}.log"
    conda:
        "../../{}/environment.yml".format(config['LDSC_dir'])
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_DBSLMM/polygenic_score_file_creator_DBSLMM.R "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input.sumstats} "
        "--plink {config[plink1_9]} "
        "--memory 5000 "
        "--ld_blocks {config[LD_block_dir]}/{wildcards[ancestry]} "
        "--rscript Rscript "
        "--dbslmm ./{config[DBSLMM_dir]}/software "
        "--munge_sumstats {input.munge_sumstats_software} "
        "--ldsc {input.ldsc_software} "
        "--ldsc_ref {config[LD_ref_dir]}/sblup_dbslmm/1000G/precomputed/{wildcards[ancestry]} "
        "--hm3_snplist {config[HapMap3_snplist_dir]}/w_hm3.snplist "
        "--sample_prev NA " # TODO
        "--pop_prev NA " # TODO
        "--output {config[Base_sumstats_dir]}/dbslmm/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[ancestry]}.{wildcards[study]} "
        "--ref_pop_scale {config[Geno_1KG_dir]}/super_pop_keep.list "
        ") &> {log}"


rule all_dbslmm_prep:
    # Run this rule to run the dbslmm method
    input: 
        expand("{}/dbslmm/{{study.study_id}}/1KGPhase3.w_hm3.{{study.ancestry}}.{{study.study_id}}.{{study.ancestry}}.scale".format(config['Base_sumstats_dir']), study=studies.itertuples())
