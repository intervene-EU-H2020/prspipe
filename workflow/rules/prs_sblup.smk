# Rules for preparing score and scale files for polygenic scoring using sblup

rule sblup_prep:
    # Implements the sblup method
    # Note that LDSC requires python 2 so Snakemake will setup this environment using the given .yml file
    # Uses the precomputed LD ref (based on 1000G) by default - 
    # If you want to compute the LD ref from scratch, replace ld_ref with expand("{}/sblup_dbslmm/1000G/fromscratch/{{study.ancestry}}/{{chr_id}}.l2.ldscore.gz".format(config['LD_ref_dir']), chr_id=range(1,23), study=studies.itertuples())
    input: 
        ld_ref=expand("{}/sblup_dbslmm/1000G/precomputed/{{study.ancestry}}/{{chr_id}}.l2.ldscore.gz".format(config['LD_ref_dir']), chr_id=range(1,23), study=studies.itertuples()),
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        super_pop_keep=rules.create_ancestry.output['super_pop_keep']
    output:
        "{}/Score_files_for_polygenic/sblup/{{study}}/1KGPhase3.w_hm3.{{study}}.{{ancestry}}.scale".format(config['Geno_1KG_dir'])
    log:
        "logs/prs_sblup_{study}.{ancestry}.log"
    conda:
        "../../{}/environment.yml".format(config['LDSC_dir'])
    threads:
        16
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_SBLUP/polygenic_score_file_creator_SBLUP.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input.qc_stats} "
        "--plink {config[plink1_9]} "
        "--gcta {config[gcta]} "
        "--munge_sumstats {config[LDSC_dir]}/munge_sumstats.py "
        "--ldsc {config[LDSC_dir]}/ldsc.py "
        "--ldsc_ref {config[LD_ref_dir]}/sblup_dbslmm/1000G/precomputed/{wildcards[ancestry]} "
        "--hm3_snplist {config[HapMap3_snplist_dir]}/w_hm3.snplist "
        "--memory 50000 "
        "--n_cores {threads} "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/sblup/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input.super_pop_keep} "
        ") &> {log}"


rule all_sblup_prep:
    # Run this rule to run the sblup method
    input: 
        expand("{}/Score_files_for_polygenic/sblup/{{study.study_id}}/1KGPhase3.w_hm3.{{study.study_id}}.{{study.ancestry}}.scale".format(config['Geno_1KG_dir']), study=studies.itertuples())
