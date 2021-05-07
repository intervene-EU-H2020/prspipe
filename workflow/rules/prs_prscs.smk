
rule download_ld_blocks:
    # 4.4.1 Prepare PRScs LD reference
    # TODO: instead of generating the reference ourselved, Oliver Pain indicated that we might want to use the pre-computed matrices instead!
    # Previously defined LD blocks have been downloaded from here: https://bitbucket.org/nygcresearch/ldetect-data/src/master/
    output:
        afr_chr=["resources/ldetect-data/AFR/fourier_ls-chr{}.bed".format(x) for x in range(1,23)],
        afr_all="resources/ldetect-data/AFR/fourier_ls-all.bed",
        asn_chr=["resources/ldetect-data/ASN/fourier_ls-chr{}.bed".format(x) for x in range(1,23)],
        ans_all="resources/ldetect-data/ASN/fourier_ls-all.bed",
        eur_chr=["resources/ldetect-data/EUR/fourier_ls-chr{}.bed".format(x) for x in range(1,23)],
        eur_all="resources/ldetect-data/EUR/fourier_ls-all.bed"
    log:
        "logs/download_ld_blocks.log"
    shell:
        "("
        "if [ -d resources/ldetect-data ]; then "
        "rm -r resources/ldetect-data ; "
        "fi;"
        "cd resources && git clone https://bitbucket.org/nygcresearch/ldetect-data.git "
        ") &> {log}"
        

rule generate_ldblock_snplist:
    # 4.4.1
    # First, we need to create a SNPlist for each LD block.
    # TODO: implement rule for all ancestries (?) (needs change in the R script)
    # TODO: the Rscript contains a hard-coded MAF and missingness filter, which should be replaced with paramters
    input:
        bed_bim_fam=rules.extract_hm3.output,
        afreq=expand(rules.allele_freq_pop.output['pop_list'], popul=['EUR'], chr=range(1,23)),
        ld_blocks=rules.download_ld_blocks.output['eur_all']
    output:
        blk_chr='{}/PRScs_LD_matrix/LD_Blocks/blk_chr'.format(config['Geno_1KG_dir']),
        blk_size='{}/PRScs_LD_matrix/LD_Blocks/blk_size'.format(config['Geno_1KG_dir']),
        snpinfo_1kg_hm3='{}/PRScs_LD_matrix/EUR/snpinfo_1kg_hm3'.format(config['Geno_1KG_dir'])
    log:
        "logs/generate_ldblock_snplist"
    shell:
        "("
        "Rscript workflow/scripts/R/prs_prscs/generate_ldblock_snplist.R "
        ") &> {log}"
        

rule calculate_ldblock:
    # 4.4.1
    # TODO: implement rule for all ancestries (?)
    input:
        rules.generate_ldblock_snplist.output
    params:
        study_ancestry="EUR"
    log:
        "logs/calculate_ldblock.log"
    output:
        touch("{}/PRScs_LD_matrix/LD_Blocks/1KGPhase3.w_hm3.EUR.all_ok".format(config['Geno_1KG_dir']))
    shell:
        "("
        "nblock=\"$(ls -1  {config[Geno_1KG_dir]}/PRScs_LD_matrix/LD_Blocks/Block_*.snplist | wc -l )\" ; "
        "echo ${{nblock}}; "
        "for block in $(seq 1 ${{nblock}} ); do"
        " {config[plink1_9]} "
        "--bfile {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--extract {config[Geno_1KG_dir]}/PRScs_LD_matrix/LD_Blocks/Block_${{block}}.snplist "
        "--r square "
        "--keep {config[Geno_1KG_dir]}/keep_files/{params[study_ancestry]}_samples.keep "
        "--out {config[Geno_1KG_dir]}/PRScs_LD_matrix/LD_Blocks/1KGPhase3.w_hm3.EUR.Block_${{block}}; done "
        ") &> {log} "
        
        
rule ldblock_to_hdf5:
    # 4.4.1
    # TODO: implement rule for all ancestries (?)
    input:
        rules.calculate_ldblock.output
    output:
        expand("{geno1kdir}/PRScs_LD_matrix/EUR/ldblk_1kg_chr{chr}.hdf5", geno1kdir=config['Geno_1KG_dir'], chr=range(1,23))
    conda:
        "../envs/py2.yaml"
    log:
        "logs/ldblock_to_hdf5.log"
    shell:
        "("
        "export HDF5_USE_FILE_LOCKING='FALSE'; "
        "python {config[GenoPred_dir]}/Scripts/PRScs_ldblk/write_ldblk_1KG_EUR.py "
        "-blkdir {config[Geno_1KG_dir]}/PRScs_LD_matrix "
        "-out {config[Geno_1KG_dir]}/PRScs_LD_matrix/EUR "
        ") &> {log}"
        
        

rule download_prscs_ld_reference:
    # Use this rule to skip the steps above and use pre-computed LD panels (for 1000 Genomes)
    # TODO: support for UKBB ? (https://github.com/getian107/PRScs)
    # TODO: maybe move to setup rules (?)
    output:
        '{}/PRScs_LD_matrix_precomputed/{{popul}}/ldblk_1kg.tar.gz'.format(config['Geno_1KG_dir']) 
    params:
        download_link=lambda wc: {
        'EUR': 'https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0',
        'AFR': 'https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=0',
        'EAS': 'https://www.dropbox.com/s/7ek4lwwf2b7f749/ldblk_1kg_eas.tar.gz?dl=0',
        }[wc['popul']]
    wildcard_constraints:
        popul='[A-Z]+'
    log:
        "logs/download_prscs_ld_reference_{popul}.log"
    shell:
        "("
        "wget -O {output} {params[download_link]}"
        ") &> {log}"
        
        
rule unpack_prscs_ld_reference_1:
    # unpack the reference (1)
    # have to wrap this with the rule below to handle different output file names
    input:
        rules.download_prscs_ld_reference.output
    output:
        snpinfo='{}/PRScs_LD_matrix_precomputed/{{popul}}/ldblk_1kg_{{popul_lower}}/snpinfo_1kg_hm3'.format(config['Geno_1KG_dir']),
        chr=expand('{geno1kdir}/PRScs_LD_matrix_precomputed/{{popul}}/ldblk_1kg_{{popul_lower}}/ldblk_1kg_chr{chr}.hdf5', geno1kdir=config['Geno_1KG_dir'], chr=range(1,23))
    log:
        "logs/unpack_prscs_ld_reference_{popul}_{popul_lower}.log"
    shell:
        "("
        "cd {config[Geno_1KG_dir]}/PRScs_LD_matrix_precomputed/{wildcards[popul]} && tar -xvf \"$(basename {input})\""
        ") &> {log}"


rule unpack_prscs_ld_reference_2:
    # unpack the reference (2)
    # example rule, not needed for workflow
    input:
        lambda wc: expand(rules.unpack_prscs_ld_reference_1.output, popul=wc['popul'], popul_lower=wc['popul'].lower())
    output:
        touch('{}/PRScs_LD_matrix_precomputed/{{popul}}/prscs_ld_reference.ok'.format(config['Geno_1KG_dir']))
        
        
rule run_prscs_precomputed_1kg:
    # 4.4.2 Run PRScs
    # runs PRScs for the 1000G populations
    # Note: unlike the original genopred, we use downloaded reference panel here
    # TODO: rule that runs it for custom reference panel
    # TODO: replace hardcoded parameters with configurable ones
    # TODO: output files
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        ld_reference=lambda wc: expand(rules.unpack_prscs_ld_reference_1.output, popul=studies.ancestry[studies.study_id == wc.study], popul_lower=studies.ancestry[studies.study_id == wc.study].str.lower())
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        ld_reference_dir = lambda wc, input: '/'.join(input['ld_reference'][0].split('/')[:-1])
    output:
        touch('{study}_PRScs_1kg.ok')
    log:
        "logs/prs_PRScs/prscs_precomputed_1kg_{study}.log"
    conda:
        "../envs/py2.yaml"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_PRScs/polygenic_score_file_creator_PRScs.R "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{params[study_ancestry]}_samples.keep "
        "--sumstats {input[qc_stats]} "
        "--plink {config[plink1_9]} "
        "--memory 5000 "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/PRScs/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--PRScs_path {config[PRScs_dir]} "
        "--PRScs_ref_path {params[ld_reference_dir]} "
        "--n_cores 8 "
        "--phi_param 1e-6,1e-4,1e-2,1,auto "
        ") &> {log} "
        
        