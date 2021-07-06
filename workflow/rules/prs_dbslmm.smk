# Rules for preparing score and scale files for polygenic scoring using dbslmm

rule download_dbslmm_ld_block:
    # Download the LD block data for the given ancestry (note that it currently only supports AFR, ASN and EUR)
    output:
        ld_block="{}/{{ancestry}}/fourier_ls-all.bed".format(config['LD_block_dir'])
    shell:
        "cd resources ;"
        "git clone https://bitbucket.org/nygcresearch/ldetect-data.git && cd ldetect-data && git checkout {config[ldbloc_version]}"


rule dbslmm_prep:
    # Implements the dbslmm method
    # Note that LDSC requires python 2 so Snakemake will setup this environment using the given .yml file
    # Uses the precomputed LD ref (based on 1000G) by default - 
    # If you want to compute the LD ref from scratch, replace ld_ref with expand("{}/sblup_dbslmm/1000G/fromscratch/{{study.ancestry}}/{{chr_id}}.l2.ldscore.gz".format(config['LD_ref_dir']), chr_id=range(1,23), study=studies.itertuples())
    # TODO: DBSLMM has a "threads" argument, which is not used in polygenic_score_file_creator_DBSLMM.R -> could potentially be used to speed things up
    input: 
        ld_ref=expand("{}/sblup_dbslmm/1000G/precomputed/{{study.ancestry}}/{{chr_id}}.l2.ldscore.gz".format(config['LD_ref_dir']), chr_id=range(1,23), study=studies.itertuples()),
        ld_block=rules.download_dbslmm_ld_block.output.ld_block,
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        super_pop_keep=rules.create_ancestry.output['super_pop_keep']
    output:
        "{}/Score_files_for_polygenic/dbslmm/{{study}}/1KGPhase3.w_hm3.{{study}}.{{ancestry}}.scale".format(config['Geno_1KG_dir'])
    log:
        "logs/prs_dbslmm_{study}.{ancestry}.log"
    conda:
        "../../{}/environment.yml".format(config['LDSC_dir'])
    threads:
        1
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_DBSLMM/polygenic_score_file_creator_DBSLMM.R "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input.qc_stats} "
        "--plink {config[plink1_9]} "
        "--memory 5000 "
        "--ld_blocks {config[LD_block_dir]}/{wildcards[ancestry]} "
        "--rscript Rscript "
        "--dbslmm ./{config[DBSLMM_dir]}/software "
        "--munge_sumstats {config[LDSC_dir]}/munge_sumstats.py "
        "--ldsc {config[LDSC_dir]}/ldsc.py "
        "--ldsc_ref {config[LD_ref_dir]}/sblup_dbslmm/1000G/precomputed/{wildcards[ancestry]} "
        "--hm3_snplist {config[HapMap3_snplist_dir]}/w_hm3.snplist "
        "--sample_prev NA " # TODO
        "--pop_prev NA " # TODO
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/dbslmm/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input.super_pop_keep} "
        ") &> {log}"


rule all_dbslmm_prep:
    # Run this rule to run the dbslmm method
    input: 
        expand("{}/Score_files_for_polygenic/dbslmm/{{study.study_id}}/1KGPhase3.w_hm3.{{study.study_id}}.{{study.ancestry}}.scale".format(config['Geno_1KG_dir']), study=studies.itertuples())
