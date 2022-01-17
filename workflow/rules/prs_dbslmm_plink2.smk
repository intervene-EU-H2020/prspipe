
#############################################################
# DBSLMM score file creator and prediction for 1000 Genomes #
#############################################################

# Rules for preparing score and scale files for polygenic scoring using dbslmm
# plink2 version should replace plink1 version.

rule download_dbslmm_ld_block:
    # Download the LD block data for the given ancestry (note that it currently only supports AFR, ASN and EUR)
    output:
        ld_block="resources/ldetect-data/{ancestry}/fourier_ls-all.bed"
    shell:
        "cd resources && "
        "git clone https://bitbucket.org/nygcresearch/ldetect-data.git"


# because the rule uses both R and ldsc + dbslmm, this makes it hard to handle... maybe create a separate image (?)

rule dbslmm_prep:
    # Implements the dbslmm method
    # Note that LDSC requires python 2 so Snakemake will setup this environment using the given .yml file
    # Uses the precomputed LD ref (based on 1000G) by default.
    # TODO: DBSLMM has a "threads" argument, which is not used in polygenic_score_file_creator_DBSLMM.R (?) -> could potentially be used to speed things up
    input:
        ld_ref=lambda wc: expand("{ldrefdir}/sblup_dbslmm/1000G/precomputed/{ancestry}/{chr}.l2.ldscore.gz", ldrefdir = config['LD_ref_dir'], ancestry=studies.ancestry[studies.study_id == wc.study].iloc[0], chr=range(1,23)),
        ld_block=lambda wc: expand(rules.download_dbslmm_ld_block.output.ld_block, ancestry=studies.ancestry[studies.study_id == wc.study].iloc[0]),
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study].iloc[0], allow_missing=True),
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        hm3_snplist=rules.download_hapmap3_snplist.output,
        geno=expand(rules.extract_hm3.output.bim, chr=range(1,23))
    output:
        touch('prs/dbslmm/{study}/ok'),
        score='prs/dbslmm/{study}/1KGPhase3.w_hm3.{study}.score.gz',
        scale=expand('prs/dbslmm/{{study}}/1KGPhase3.w_hm3.{{study}}.{ancestry}.scale', ancestry=config['1kg_superpop'])
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        ld_ref_dir=lambda wc, input: '/'.join(input['ld_ref'][0].split('/')[-1])
    log:
        "logs/prs_dbslmm_{study}.log"
    conda:
        "../envs/ldsc.yaml"
    threads:
        1
    resources:
        mem_mb=5000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_DBSLMM/polygenic_score_file_creator_DBSLMM_plink2.R "
        "--ref_plink_chr resources/1kg/1KGPhase3.w_hm3.chr "
        "--ref_keep resources/1kg/keep_files/{params[study_ancestry]}_samples.keep "
        "--sumstats {input.qc_stats} "
        "--plink {config[plink1_9]} "
        "--memory {resources[mem_mb]} "
        "--ld_blocks resources/ldetect-data/{params[study_ancestry]}/ "
        "--rscript Rscript "
        "--dbslmm ./{config[DBSLMM_dir]}/software "
        "--munge_sumstats {config[LDSC_dir]}/munge_sumstats.py "
        "--ldsc {config[LDSC_dir]}/ldsc.py "
        "--ldsc_ref {config[LD_ref_dir]}/sblup_dbslmm/1000G/precomputed/{params[study_ancestry]} "
        "--hm3_snplist {input.hm3_snplist} "
        "--sample_prev NA " # TODO
        "--pop_prev NA " # TODO
        "--output prs/dbslmm/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input.super_pop_keep} "
        "--ignore OR && " # TODO -> should this always be OR, or can it also be BETA?
        "gzip prs/dbslmm/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.score "
        ") &> {log}"

rule all_dbslmm_prep:
    # Run this rule to run the dbslmm method
    input: 
        expand(rules.dbslmm_prep.output, study=studies.study_id)
