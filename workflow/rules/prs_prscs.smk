
#TODO: change LD matrix output directories, they're a mess... (also for other methods)

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
        blk_chr='resources/LD_matrix/prscs/1000G/fromscratch/EUR/LD_Blocks/blk_chr',
        blk_size='resources/LD_matrix/prscs/1000G/fromscratch/EUR/LD_Blocks/blk_size',
        snpinfo_1kg_hm3='resources/LD_matrix/prscs/1000G/fromscratch/EUR/snpinfo_1kg_hm3'
    log:
        "logs/generate_ldblock_snplist.log"
    shell:
        "("
        "Rscript workflow/scripts/R/prs_prscs/generate_ldblock_snplist.R EUR"
        ") &> {log}"
        

rule calculate_ldblock:
    # 4.4.1
    # TODO: implement rule for all ancestries (?)
    # TODO: output files!
    input:
        rules.generate_ldblock_snplist.output
    params:
        study_ancestry="EUR"
    log:
        "logs/calculate_ldblock.log"
    output:
        # creates many files with '.ld' suffix...
        touch('resources/LD_matrix/prscs/1000G/fromscratch/EUR/LD_Blocks/ldblk.ok')
    shell:
        "("
        "nblock=\"$(ls -1  ./resources/LD_matrix/prscs/1000G/fromscratch/EUR/LD_Blocks/Block_*.snplist | wc -l )\" ; "
        "echo ${{nblock}}; "
        "for block in $(seq 1 ${{nblock}} ); do "
        "{config[plink1_9]} "
        "--bfile {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--extract resources/LD_matrix/prscs/1000G/fromscratch/EUR/LD_Blocks/Block_${{block}}.snplist "
        "--r square "
        "--keep {config[Geno_1KG_dir]}/keep_files/{params[study_ancestry]}_samples.keep "
        "--out resources/LD_matrix/prscs/1000G/fromscratch/EUR/LD_blocks/1KGPhase3.w_hm3.EUR.Block_${{block}}; done "
        ") &> {log} "
        
        
rule ldblock_to_hdf5:
    # 4.4.1
    # TODO: implement rule for all ancestries (?)
    input:
        rules.calculate_ldblock.output
    output:
        expand("resources/LD_matrix/prscs/1000G/fromscratch/EUR/ldblk_1kg_chr{chr}.hdf5", chr=range(1,23))
    conda:
        "../envs/py2.yaml"
    log:
        "logs/ldblock_to_hdf5.log"
    shell:
        "("
        "export HDF5_USE_FILE_LOCKING='FALSE'; "
        "python {config[GenoPred_dir]}/Scripts/PRScs_ldblk/write_ldblk_1KG_EUR.py "
        "-blkdir resources/LD_matrix/prscs/1000G/fromscratch/EUR "
        "-out resources/LD_matrix/prscs/1000G/fromscratch/EUR "
        ") &> {log}"
        
        

rule download_prscs_ld_reference:
    # Use this rule to skip the steps above and use pre-computed LD panels (for 1000 Genomes)
    output:
        'resources/LD_matrix/prscs/1000G/precomputed/{popul}/ldblk_1kg.tar.gz' 
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


rule download_prscs_ld_reference_ukbb:
    # Use this rule to skip the steps above and use pre-computed LD panels (for UK Biobank)
    output:
        'resources/LD_matrix/prscs/UKBB/precomputed/{popul}/ldblk_ukbb.tar.gz' 
    params:
        download_link=lambda wc: {
        'EUR': 'https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=0',
        'AFR': 'https://www.dropbox.com/s/dtccsidwlb6pbtv/ldblk_ukbb_afr.tar.gz?dl=0',
        'EAS': 'https://www.dropbox.com/s/fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz?dl=0',
        'SAS': 'https://www.dropbox.com/s/nto6gdajq8qfhh0/ldblk_ukbb_sas.tar.gz?dl=0'
        }[wc['popul']]
    wildcard_constraints:
        popul='[A-Z]+'
    log:
        "logs/download_prscs_ld_reference_ukbb_{popul}.log"
    shell:
        "("
        "wget -O {output} {params[download_link]}"
        ") &> {log}"
        
        
rule unpack_prscs_ld_reference_1:
    # unpack the reference for 1000 Genomes
    # note: the path has to contain "1kg" or "ukbb", otherwise PRScs crashes
    input:
        rules.download_prscs_ld_reference.output
    output:
        snpinfo='resources/LD_matrix/prscs/1000G/precomputed/{popul}/ldblk_1kg/snpinfo_1kg_hm3',
        chr=expand('resources/LD_matrix/prscs/1000G/precomputed/{popul}/ldblk_1kg/ldblk_1kg_chr{chr}.hdf5', chr=range(1,23), allow_missing=True)
    log:
        "logs/unpack_prscs_ld_reference_{popul}.log"
    shell:
        "("
        "cd resources/LD_matrix/prscs/1000G/precomputed/{wildcards[popul]} && tar -xvf \"$(basename {input})\" && "
        "mv ldblk_1kg_*/ ./ldblk_1kg && rm \"$(basename {input})\" "
        ") &> {log}"
        
        
rule unpack_prscs_ld_reference_1_ukbb:
    # unpack the reference for UK Biobank
    # note: the path has to contain "1kg" or "ukbb", otherwise PRScs crashes
    input:
        rules.download_prscs_ld_reference_ukbb.output
    output:
        snpinfo='resources/LD_matrix/prscs/UKBB/precomputed/{popul}/ldblk_ukbb/snpinfo_ukbb_hm3',
        chr=expand('resources/LD_matrix/prscs/UKBB/precomputed/{popul}/ldblk_ukbb/ldblk_ukbb_chr{chr}.hdf5', chr=range(1,23), allow_missing=True)
    log:
        "logs/unpack_prscs_ld_reference_ukbb_{popul}.log"
    shell:
        "("
        "cd resources/LD_matrix/prscs/UKBB/precomputed/{wildcards[popul]} && tar -xvf \"$(basename {input})\" && "
        "mv ldblk_ukbb_*/ ./ldblk_ukbb && rm \"$(basename {input})\" "
        ") &> {log}"


#TODO: rule that runs PRScs with the custom LD reference panel


rule run_prscs_precompld_1kg:
    # 4.4.2 Run PRScs with precomputed 1000G LD reference panel
    # runs PRScs for the 1000G populations
    # TODO: replace hardcoded parameters with configurable ones
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        ld_reference=lambda wc: expand(rules.unpack_prscs_ld_reference_1.output, popul=studies.ancestry[studies.study_id == wc.study])
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        ld_reference_dir = lambda wc, input: '/'.join(input['ld_reference'][0].split('/')[:-1])
    output:
        adj_ss=expand('{geno1kg}/Score_files_for_polygenic/PRScs_precompld/{{study}}/1KGPhase3.w_hm3.{{study}}_pst_eff_a1_b0.5_phi{phi}_chr8.txt', geno1kg=config['Geno_1KG_dir'], phi=['1e+00','1e-02','1e-04','1e-06','auto'], chr=range(1,23)),
        scale=expand('{geno1kg}/Score_files_for_polygenic/PRScs_precompld/{{study}}/1KGPhase3.w_hm3.{{study}}.{superpop}.scale', geno1kg=config['Geno_1KG_dir'], superpop=config['1kg_superpop']),
        log='{}/Score_files_for_polygenic/PRScs_precompld/{{study}}/1KGPhase3.w_hm3.{{study}}.log'.format(config['Geno_1KG_dir']),
        # touch('{}/Score_files_for_polygenic/PRScs_precompld/{{study}}/prscs.ok'.format(config['Geno_1KG_dir']))
    log:
        "logs/prs_PRScs/prscs_precompld_1kg_{study}.log"
    conda:
        "../envs/py2.yaml"
    threads:
        16
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_PRScs/polygenic_score_file_creator_PRScs.R "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{params[study_ancestry]}_samples.keep "
        "--sumstats {input[qc_stats]} "
        "--plink {config[plink1_9]} "
        "--memory 5000 "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/PRScs_precompld/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--PRScs_path {config[PRScs_dir]} "
        "--PRScs_ref_path {params[ld_reference_dir]} "
        "--n_cores {threads} "
        "--phi_param 1e-6,1e-4,1e-2,1,auto "
        ") &> {log} "
        

rule run_prscs_precompld_1kg_refukbb:
    # 4.4.2 Run PRScs with precomputed UKBB LD reference panel
    # runs PRScs for the 1000G populations
    # TODO: replace hardcoded parameters with configurable ones
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        ld_reference=lambda wc: expand(rules.unpack_prscs_ld_reference_1_ukbb.output, popul=studies.ancestry[studies.study_id == wc.study])
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        ld_reference_dir = lambda wc, input: '/'.join(input['ld_reference'][0].split('/')[:-1])
    output:
        adj_ss=expand('{geno1kg}/Score_files_for_polygenic/PRScs_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}_pst_eff_a1_b0.5_phi{phi}_chr8.txt', geno1kg=config['Geno_1KG_dir'], phi=['1e+00','1e-02','1e-04','1e-06','auto'], chr=range(1,23)),
        scale=expand('{geno1kg}/Score_files_for_polygenic/PRScs_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.{superpop}.scale', geno1kg=config['Geno_1KG_dir'], superpop=config['1kg_superpop']),
        log='{}/Score_files_for_polygenic/PRScs_precompld_ukbb/{{study}}/1KGPhase3.w_hm3.{{study}}.log'.format(config['Geno_1KG_dir'])
    log:
        "logs/prs_PRScs/prscs_precompld_1kg_refukbb_{study}.log"
    conda:
        "../envs/py2.yaml"
    threads:
        32
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_PRScs/polygenic_score_file_creator_PRScs.R "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{params[study_ancestry]}_samples.keep "
        "--sumstats {input[qc_stats]} "
        "--plink {config[plink1_9]} "
        "--memory 5000 "
        "--output {config[Geno_1KG_dir]}/Score_files_for_polygenic/PRScs_precompld_ukbb/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--PRScs_path {config[PRScs_dir]} "
        "--PRScs_ref_path {params[ld_reference_dir]} "
        "--n_cores {threads} "
        "--phi_param 1e-6,1e-4,1e-2,1,auto "
        ") &> {log} "