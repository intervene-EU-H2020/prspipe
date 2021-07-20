
rule download_genetic_map:
    # 4.6 Prepare score and scale files for polygenic scoring using SBayesR
    # somehow these downloads take forever...
    output:
        expand("{genetic_map}/CEU/chr{chr}.OMNI.interpolated_genetic_map", genetic_map=config['genetic_map'], chr=range(1,23))
    log:
        "logs/download_genetic_map.log"
    shell:
        "("
        "cd {config[genetic_map]}/CEU; "
        "for chr in $(seq 1 22); do "
        "wget https://github.com/joepickrell/1000-genomes-genetic-maps/raw/master/interpolated_OMNI/chr${{chr}}.OMNI.interpolated_genetic_map.gz; "
        "done; "
        "gunzip *.gz "
        ") &> {log} "
        

rule generate_ldm_chunks_1:
    # generate sbayesr LD matrix in chunks, step 1
    # filter SNPs to those with MAF > 0.0
    # generate chunk files
    input:
        rules.run_allele_freq_superpop.output
    output:
        # A rule with dynamic output may not define any non-dynamic output files.
        # however, this temporary file is produced as well:
        # snp=expand("resources/LD_matrix/sbayesr/1000G/fromscratch/{{popul}}/chr{{chr}}/SNP.txt", geno1kdir=config['Geno_1KG_dir']),
        chunk=dynamic(expand("resources/LD_matrix/sbayesr/1000G/fromscratch/{{popul}}/chr{{chr}}/SNP_{{chunk}}.txt", geno1kdir=config['Geno_1KG_dir']))
    shell:
        "Rscript workflow/scripts/R/prs_sbayesr/generate_ld_matrix_chunks.R {wildcards[popul]} {config[sbayesr_ldm_chunksize]} {wildcards[chr]}"


rule generate_ldm_chunks_1_all:
    # example of how to run the dynamic rule above
    input:
        dynamic(expand(rules.generate_ldm_chunks_1.output, popul='EUR', chr=range(1,23), allow_missing=True))


rule generate_ldm_chunks_2:
    # generate sbayesr LD matrix in chunks, step 2
    input:
        chunkfile='resources/LD_matrix/sbayesr/1000G/fromscratch/{popul}/chr{chr}/SNP_{chunk}.txt',
        genetic_map=lambda wc: rules.download_genetic_map.output[int(wc['chr'])-1],
        extract_hm3=rules.extract_hm3.output,
        keep=rules.create_ancestry.output
    output:
        bin=temp('resources/LD_matrix/sbayesr/1000G/fromscratch/{popul}/chr{chr}/1KGPhase3.w_hm3.{chunk}.ldm.shrunk.bin'),
        info='resources/LD_matrix/sbayesr/1000G/fromscratch/{popul}/chr{chr}/1KGPhase3.w_hm3.{chunk}.ldm.shrunk.info'
    log:
        'logs/generate_ldm_chunks_2/{popul}_{chr}_{chunk}.log'
    shell:
        "("
        "while read start end; do "
        "{config[gctb]} "
        "--bfile {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr{wildcards[chr]} "
        "--keep {config[Geno_1KG_dir]}/keep_files/{wildcards[popul]}_samples.keep "
        "--make-shrunk-ldm "
        "--extract resources/LD_matrix/sbayesr/1000G/fromscratch/{wildcards[popul]}/chr{wildcards[chr]}/SNP.txt "
        "--gen-map {input[genetic_map]} "
        "--snp ${{start}}-${{end}} "
        "--out resources/LD_matrix/sbayesr/1000G/fromscratch/{wildcards[popul]}/chr{wildcards[chr]}/1KGPhase3.w_hm3; "
        "done < {input[chunkfile]} "
        ") &> {log} "


rule merge_ldm_chunks:
    # merge the chunks generated by the rule above
    input:
        chunks=dynamic(rules.generate_ldm_chunks_2.output.bin)
    output:
        'resources/LD_matrix/sbayesr/1000G/fromscratch/{popul}/1KGPhase3.w_hm3.{popul}.chr{chr}.ldm.shrunk.bin'
    log:
        "logs/merge_ldm_chunks/{popul}_{chr}.log"
    shell:
        "("
        "Rscript workflow/scripts/R/prs_sbayesr/merge_ld_matrix_chunks.R {wildcards[popul]} {wildcards[chr]} {input[chunks]} "
        ") &> {log} "
        

rule all_sbayesr_ldm:
    # run the rule above for all chromosomes and EUR reference population
    input:
        expand(rules.merge_ldm_chunks.output, popul=['EUR'], chr=range(1,23))
        

rule download_sbasesr_ld_reference:
    # download the pre-computed LD matrix provided for GCTB (EUR UKBB)
    # effectively this skips the steps to generate the LD-matrix (ldm) above
    # TODO: currently only works for EUR, and path is probably hardcoded in subsequent rules
    # Genopred 4.6
    output:
        bin=expand('resources/LD_matrix/sbayesr/UKBB/precomputed/EUR/ukb{popul}u_hm3_v3_50k_chr{chr}.ldm.sparse.bin', geno1kdir=config['Geno_1KG_dir'], chr=range(1,23), allow_missing=True),
        info=expand('resources/LD_matrix/sbayesr/UKBB/precomputed/EUR/ukb{popul}u_hm3_v3_50k_chr{chr}.ldm.sparse.info', geno1kdir=config['Geno_1KG_dir'],chr=range(1,23), allow_missing=True)
    log:
        'logs/download_sbasesr_ld_reference_{popul}.log'
    shell:
        "("
        "bash workflow/scripts/bash/download_sbayesr_reference.sh {wildcards[popul]} resources/LD_matrix/sbayesr/1000G/precomputed/EUR/ "
        ") &> {log}"
        
  
rule all_download_sbayesr_ld_reference:
    # runs rule above for all supported superpopulations (currently only EUR)
    input:
        expand(rules.download_sbasesr_ld_reference.output, popul=['EUR'])
        
        
# TODO: rule that runs SBayesR with custom reference panel (1000 genomes)
        
rule run_sbayesr_precompld_1kg_refukbb_robust:
    # 4.6 Prepare score and scale files for polygenic scoring using SBayesR
    # rule that uses pre-computed LD matrix
    # this uses the --robust setting
    # TODO: currently assumes EUR ancestry -> could change in the future
    # TODO: implement the rule for LD-matrix that is computed "from scratch"
    # TODO: get rid of ugly names (?)
    # TODO: adjust outputs
    # TODO: this only runs SBayesR with default parameters (?, except --rsq 0.95 ), what about other parameters
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        ld_reference=lambda wc: expand(rules.download_sbasesr_ld_reference.output, popul=studies.ancestry[studies.study_id == wc.study])
    output:
        log1=expand('{geno1kg}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{{study}}/SBayesR{setting}.chr{chr}.log', geno1kg=config['Geno_1KG_dir'], setting=['','.robust'], chr=range(1,23)),
        log2=expand('{geno1kg}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.log', geno1kg=config['Geno_1KG_dir']),
        scale=expand('{geno1kg}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.{superpop}.scale', geno1kg=config['Geno_1KG_dir'], superpop=config['1kg_superpop']),
        parres=expand('{geno1kg}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{{study}}/GWAS_sumstats_SBayesR{setting}.GW.parRes', geno1kg=config['Geno_1KG_dir'], setting=['','.robust']),
        snpres=expand('{geno1kg}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{{study}}/GWAS_sumstats_SBayesR{setting}.GW.snpRes', geno1kg=config['Geno_1KG_dir'], setting=['','.robust'])
    log:
        "logs/run_sbayesr_precomputed_1kg_{study}.log"
    threads:
        22
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_SBayesR/polygenic_score_file_creator_SBayesR.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--sumstats {input[qc_stats]} "
        "--plink {config[plink1_9]} "
        "--gctb {config[gctb]} "
        "--ld_matrix_chr resources/LD_matrix/sbayesr/UKBB/precomputed/EUR/ukbEURu_hm3_v3_50k_chr "
        "--memory 50000 "
        "--n_cores {threads} "
        "--robust TRUE "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {config[Geno_1KG_dir]}/super_pop_keep.list "
        ") &> {log}"
        

rule all_run_sbayesr_precompld_1kg_refukbb_robust:
    # runs rule above for all studies
    input:
        expand(rules.run_sbayesr_precompld_1kg_refukbb_robust.output, study=studies.study_id)
        

# robust analysis and "regular" analysis are performed in the same rule
# the rules below are not needed.
#rule run_sbayesr_precompld_1kg_refukbb:
#    # 4.6 Prepare score and scale files for polygenic scoring using SBayesR
#    # rule that uses pre-computed LD matrix
#    # this does not use the robust setting, but instead forces UKBB allele coding!
#    # TODO: currently assumes EUR ancestry -> could change in the future
#    # TODO: implement the rule for LD-matrix that is computed "from scratch"
#    # TODO: get rid of ugly names (?)
#    # TODO: adjust outputs
#    # TODO: this only runs SBayesR with default parameters (?), what about other parameters
#    input:
#        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
#        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
#        ld_reference=lambda wc: expand(rules.download_sbasesr_ld_reference.output, popul=studies.ancestry[studies.study_id == wc.study]),
#        ref_maf = lambda wc: expand(config['UKBB_output'] + '/UKBB_ref/genotype/UKBB.noPheno.{ancestry}.10K.chr{chr}.frq', chr=range(1,23), ancestry=studies.ancestry[studies.study_id == wc.study])
#        #ref_maf=rules.extract_ukbb10k.output['frq']
#    output:
#        log1=expand('{geno1kg}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{{study}}/SBayesR.chr{chr}.log', geno1kg=config['Geno_1KG_dir'], chr=range(1,23)),
#        log2=expand('{geno1kg}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.log', geno1kg=config['Geno_1KG_dir']),
#        scale=expand('{geno1kg}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.{superpop}.scale', geno1kg=config['Geno_1KG_dir'], superpop=config['1kg_superpop']),
#        parres=expand('{geno1kg}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{{study}}/GWAS_sumstats_SBayesR.GW.parRes', geno1kg=config['Geno_1KG_dir']),
#        snpres=expand('{geno1kg}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{{study}}/GWAS_sumstats_SBayesR.GW.snpRes', geno1kg=config['Geno_1KG_dir'])
#    params:
#        ancestry = lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0]
#    log:
#        "logs/run_sbayesr_precomputed_1kg_{study}.log"
#    threads:
#        22
#    shell:
#        # this argument actually doesn't do anything (?)
#        # "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[popul]}_samples.keep "
#        "("
#        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_SBayesR/polygenic_score_file_creator_SBayesR.R "
#        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
#        "--sumstats {input[qc_stats]} "
#        "--plink {config[plink1_9]} "
#        "--gctb {config[gctb]} "
#        "--ld_matrix_chr resources/LD_matrix/sbayesr/UKBB/precomputed/EUR/ukbEURu_hm3_v3_50k_chr "
#        "--memory 50000 "
#        "--n_cores {threads} "
#        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/SBayesR_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
#        "--ref_pop_scale {config[Geno_1KG_dir]}/super_pop_keep.list "
#        ") &> {log}"
#        # "--force_ref_frq results/UKBB/UKBB_ref/genotype/UKBB.noPheno.{params[ancestry]}.10K "
#        
#
#rule all_run_sbayesr_precompld_1kg_refukbb:
#    # runs rule above for all studies
#    input:
#        expand(rules.run_sbayesr_precompld_1kg_refukbb.output, study=studies.study_id)
