
rule install_gctb:
    # install GCTB to the default location
    output:
        "bin/gctb/gctb"
    shell:
        "cd bin && "
        "wget https://cnsgenomics.com/software/gctb/download/gctb_2.03beta_Linux.zip && "
        "unzip gctb_2.03beta_Linux.zip && "
        "rm gctb_2.03beta_Linux.zip && "
        "mv gctb_2.03beta_Linux/* ./gctb"


rule download_sbasesr_ld_reference:
    # download the pre-computed LD matrix provided for GCTB (EUR UKBB)
    # TODO: currently only works for EUR, and path is probably hardcoded in subsequent rules
    # Genopred 4.6
    output:
        bin=expand('resources/LD_matrix/sbayesr/UKBB/precomputed/EUR/ukb{ancestry}u_hm3_v3_50k_chr{chr}.ldm.sparse.bin', chr=range(1,23), allow_missing=True),
        info=expand('resources/LD_matrix/sbayesr/UKBB/precomputed/EUR/ukb{ancestry}u_hm3_v3_50k_chr{chr}.ldm.sparse.info',chr=range(1,23), allow_missing=True)
    log:
        'logs/download_sbasesr_ld_reference/{ancestry}.log'
    shell:
        "("
        "bash workflow/scripts/bash/download_sbayesr_reference.sh {wildcards[ancestry]} resources/LD_matrix/sbayesr/1000G/precomputed/EUR/ "
        ") &> {log}"
        
  
rule all_download_sbayesr_ld_reference:
    # runs rule above for all supported superpopulations (currently only EUR)
    input:
        expand(rules.download_sbasesr_ld_reference.output, ancestry=['EUR'])
        
                
rule prs_scoring_sbayesr:
    # 4.6 Prepare score and scale files for polygenic scoring using SBayesR and pre-computed LD matrix
    # TODO: currently assumes EUR ancestry -> could change in the future
    # TODO: this only runs SBayesR with default parameters (?, except --rsq 0.95 ), what about other parameters
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        ld_reference=lambda wc: expand(rules.download_sbasesr_ld_reference.output, ancestry=studies.ancestry[studies.study_id == wc.study]),
        hm3_gw=rules.extract_hm3_gw.output,
        gctb=config['gctb']
    output:
        touch('prs/sbayesr/{study}/ok')
    log:
        "logs/prs_scoring_sbayesr/{study}.log"
    threads:
        16
    resources:
        mem_mb=32000,
        time="04:00:00",
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_SBayesR/polygenic_score_file_creator_SBayesR_plink2.R "
        "--ref_plink resources/1kg/1KGPhase3.w_hm3.GW "
        "--sumstats {input[qc_stats]} "
        "--plink {config[plink1_9]} "
        "--gctb {config[gctb]} "
        "--ld_matrix_chr resources/LD_matrix/sbayesr/UKBB/precomputed/EUR/ukbEURu_hm3_v3_50k_chr "
        "--memory {resources[mem_mb]} "
        "--n_cores {threads} "
        "--robust T "
        "--output prs/sbayesr/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input[super_pop_keep]} "
        ") &> {log}"
        

rule all_run_sbayesr_precompld_1kg_refukbb_robust:
    # runs rule above for all studies
    input:
        expand(rules.prs_scoring_sbayesr.output, study=studies.study_id)
        
