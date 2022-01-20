# Rules for preparing score and scale files for polygenic scoring using sblup

rule install_gcta:
    # install GCTA to the default location
    output:
        "bin/gcta/gcta64"
    shell:
        "cd bin && "
        "wget https://cnsgenomics.com/software/gcta/bin/gcta_1.93.2beta.zip && "
        "unzip gcta_1.93.2beta.zip && "
        "rm gcta_1.93.2beta.zip && "
        "mv gcta_1.93.2beta/* ./gcta"

rule prs_scoring_sblup:
    # Implements the sblup method
/{{chr_id}}.l2.ldscore.gz".format(config['LD_ref_dir']), chr_id=range(1,23), study=studies.itertuples())
    input: 
        ld_ref=expand("resources/LD_matrix/sblup_dbslmm/1000G/precomputed/{{ancestry}}/{chr}.l2.ldscore.gz", chr=range(1,23)),
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        gcta=config['gcta']
    output:
        "prs/sblup/{study}/1KGPhase3.w_hm3.{study}.{ancestry}.scale"
    log:
        "logs/prs_scoring_sblup/{study}_{ancestry}.log"
    conda:
        "../envs/ldsc.yaml"
    threads:
        16
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_SBLUP/polygenic_score_file_creator_SBLUP.R "
        "--ref_plink resources/1kg/1KGPhase3.w_hm3.GW "
        "--ref_keep resources/1kg/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input.qc_stats} "
        "--plink {config[plink1_9]} "
        "--gcta {config[gcta]} "
        "--munge_sumstats {config[LDSC_dir]}/munge_sumstats.py "
        "--ldsc {config[LDSC_dir]}/ldsc.py "
        "--ldsc_ref {config[LD_ref_dir]}/sblup_dbslmm/1000G/precomputed/{wildcards[ancestry]} "
        "--hm3_snplist resources/HapMap3_snplist/w_hm3.snplist "
        "--memory 50000 "
        "--n_cores {threads} "
        "--output prs/sblup/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input.super_pop_keep} "
        ") &> {log}"


rule all_prs_scoring_sblup:
    # Run this rule to run the sblup method
    input: 
        expand(rules.prs_scoring_sblup.output, zip, study=studies.study_id, ancestry=studies.ancestry)
