
#############################################################
# DBSLMM score file creator and prediction for 1000 Genomes #
#############################################################

# Rules for preparing score and scale files for polygenic scoring using dbslmm

rule install_dbslmm:
    # TODO: define specific version
    output:
        directory('workflow/scripts/DBSLMM')
    shell:
        "mkdir -p workflow/scripts/DBSLMM && "
        "git clone https://github.com/intervene-EU-H2020/DBSLMM.git workflow/scripts/DBSLMM"
    

rule download_dbslmm_ld_block:
    # Download the LD block data for the given ancestry (note that it currently only supports AFR, ASN and EUR)
    output:
        ld_block=expand("resources/ldetect-data/{ancestry}/fourier_ls-all.bed", ancestry=['AFR','ASN','EUR'])
    shell:
        "cd resources && "
        "git clone https://bitbucket.org/nygcresearch/ldetect-data.git"


rule prs_scoring_dbslmm:
    # Implements the dbslmm method
    # Uses the precomputed LD ref (based on 1000G) by default.
    # TODO: DBSLMM has a "threads" argument, which is not used in polygenic_score_file_creator_DBSLMM.R (?) -> could potentially be used to speed things up
    input:
        ld_ref=lambda wc: expand("resources/LD_matrix/sblup_dbslmm/1000G/precomputed/{ancestry}/{chr}.l2.ldscore.gz", ancestry=studies.ancestry[studies.study_id == wc.study].iloc[0], chr=range(1,23)),
        ld_block=rules.download_dbslmm_ld_block.output,
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study].iloc[0], allow_missing=True),
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        hm3_snplist=rules.download_hapmap3_snplist.output,
        geno=expand(rules.extract_hm3.output.bim, chr=range(1,23)),
        dbslmm=rules.install_dbslmm.output
    output:
        touch('prs/dbslmm/{study}/ok'),
        score='prs/dbslmm/{study}/1KGPhase3.w_hm3.{study}.score.gz',
        scale=expand('prs/dbslmm/{{study}}/1KGPhase3.w_hm3.{{study}}.{ancestry}.scale', ancestry=config['1kg_superpop'])
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        ld_ref_dir=lambda wc, input: '/'.join(input['ld_ref'][0].split('/')[-1])
    log:
        "logs/prs_scoring_dbslmm/{study}.log"
    conda:
        "../envs/ldsc.yaml"
    threads:
        1
    resources:
        mem_mb=5000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_3.sqsh --no-container-mount-home"
    shell:
        "("
        "{config[Rscript]} {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_DBSLMM/polygenic_score_file_creator_DBSLMM_plink2.R "
        "--ref_plink_chr resources/1kg/1KGPhase3.w_hm3.chr "
        "--ref_keep resources/1kg/keep_files/{params[study_ancestry]}_samples.keep "
        "--sumstats {input.qc_stats} "
        "--plink {config[plink1_9]} "
        "--memory {resources[mem_mb]} "
        "--ld_blocks resources/ldetect-data/{params[study_ancestry]}/ "
        "--rscript {config[Rscript]} "
        "--dbslmm workflow/scripts/DBSLMM/software "
        "--munge_sumstats {config[LDSC_dir]}/munge_sumstats.py "
        "--ldsc {config[LDSC_dir]}/ldsc.py "
        "--ldsc_ref resources/LD_matrix/sblup_dbslmm/1000G/precomputed/{params[study_ancestry]} "
        "--hm3_snplist {input.hm3_snplist} "
        "--sample_prev NA " # TODO
        "--pop_prev NA " # TODO
        "--output prs/dbslmm/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input.super_pop_keep} "
        "--ignore OR && " # TODO -> should this always be OR, or can it also be BETA?
        "gzip prs/dbslmm/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]}.score "
        ") &> {log}"


rule all_prs_scoring_dbslmm:
    # Run this rule to run the dbslmm method
    input: 
        expand(rules.prs_scoring_dbslmm.output, study=studies.study_id)
