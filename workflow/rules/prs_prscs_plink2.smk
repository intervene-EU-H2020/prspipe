
rule install_prscs:
    output:
        "workflow/scripts/PRScs/PRScs.py"
    params:
        PRSCS_VERSION="f2f2b4201ffe80715d4bc46582a5207f6a31dd57"
    shell:
        "if [ -d workflow/scripts/PRScs ]; then rm -rf workflow/scripts/PRScs; fi; "
        "mkdir -p workflow/scripts/PRScs && "
        "git clone https://github.com/getian107/PRScs.git workflow/scripts/PRScs && "
        "cd workflow/scripts/PRScs && "
        "git checkout {params[PRSCS_VERSION]} "
        
    

rule download_prscs_ld_reference_ukbb:
    # Use this rule to skip the steps above and use pre-computed LD panels (for UK Biobank)
    output:
        'resources/LD_matrix/prscs/UKBB/precomputed/{ancestry}/ldblk_ukbb.tar.gz' 
    params:
        download_link=lambda wc: {
        'EUR': 'https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=0',
        'AFR': 'https://www.dropbox.com/s/dtccsidwlb6pbtv/ldblk_ukbb_afr.tar.gz?dl=0',
        'EAS': 'https://www.dropbox.com/s/fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz?dl=0',
        'SAS': 'https://www.dropbox.com/s/nto6gdajq8qfhh0/ldblk_ukbb_sas.tar.gz?dl=0'
        }[wc['ancestry']]
    wildcard_constraints:
        ancestry='[A-Z]+'
    log:
        "logs/download_prscs_ld_reference_ukbb/{ancestry}.log"
    shell:
        "("
        "wget -O {output} {params[download_link]}"
        ") &> {log}"
        
        
rule unpack_prscs_ld_reference_ukbb:
    # unpack the reference for UK Biobank
    # note: the path has to contain "1kg" or "ukbb", otherwise PRScs crashes
    input:
        rules.download_prscs_ld_reference_ukbb.output
    output:
        snpinfo='resources/LD_matrix/prscs/UKBB/precomputed/{ancestry}/ldblk_ukbb/snpinfo_ukbb_hm3',
        chr=expand('resources/LD_matrix/prscs/UKBB/precomputed/{ancestry}/ldblk_ukbb/ldblk_ukbb_chr{chr}.hdf5', chr=range(1,23), allow_missing=True)
    log:
        "logs/unpack_prscs_ld_reference_ukbb/{ancestry}.log"
    shell:
        "("
        "cd resources/LD_matrix/prscs/UKBB/precomputed/{wildcards[ancestry]} && tar -xvf \"$(basename {input})\" && "
        "if [ -d ldblk_ukbb ]; then rm -r ldblk_ukbb; fi; "
        "mv ldblk_ukbb_*/ ./ldblk_ukbb && rm \"$(basename {input})\" "        
        ") &> {log}"


rule prs_scoring_prscs:
    # 4.4.2 Run PRScs with precomputed UKBB LD reference panel
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        ld_reference=lambda wc: expand(rules.unpack_prscs_ld_reference_ukbb.output, ancestry=studies.ancestry[studies.study_id == wc.study]),
        hm3=expand(rules.extract_hm3.output, chr=range(1,23)),
        prscs=config['PRScs_dir']
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
        ld_reference_dir = lambda wc, input: '/'.join(input['ld_reference'][0].split('/')[:-1])
    output:
        touch('prs/prscs/{study}/ok'),
        score='prs/prscs/{study}/1KGPhase3.w_hm3.{study}.score.gz',
        scale=expand('prs/prscs/{{study}}/1KGPhase3.w_hm3.{{study}}.{ancestry}.scale', ancestry=config['1kg_superpop']),
        log='prs/prscs/{study}/1KGPhase3.w_hm3.{study}.log'
    log:
        "logs/prs_scoring_prscs/{study}.log"
    conda:
        "../envs/ldsc.yaml"
    threads:
        16
    resources:
        mem_mb=32000,
        time="12:00:00",
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh --no-container-mount-home"
    shell:
        "( "
        "export MKL_NUM_THREADS=1; "
        "export NUMEXPR_NUM_THREADS=1; "
        "export OMP_NUM_THREADS=1; "
        "export OPENBLAS_NUM_THREADS=1; "
        "{config[Rscript]} {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_PRScs/polygenic_score_file_creator_PRScs_plink2.R "
        "--ref_plink_chr resources/1kg/1KGPhase3.w_hm3.chr "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--sumstats {input[qc_stats]} "
        "--plink2 {config[plink2]} "
        "--memory 5000 "
        "--output prs/prscs/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--PRScs_path {config[PRScs_dir]} "
        "--PRScs_ref_path {params[ld_reference_dir]} "
        "--n_cores {threads} "
        "--phi_param 1e-6,1e-4,1e-2,1,auto "
        "--seed 1 "
        ") &> {log} "


rule all_prs_scoring_prscs:
    input:
        expand(rules.prs_scoring_prscs.output, study=studies.study_id)



# rule download_prscs_ld_reference:
#     # Use this rule to skip the steps above and use pre-computed LD panels (for 1000 Genomes)
#     output:
#         'resources/LD_matrix/prscs/1000G/precomputed/{popul}/ldblk_1kg.tar.gz' 
#     params:
#         download_link=lambda wc: {
#         'EUR': 'https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0',
#         'AFR': 'https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=0',
#         'EAS': 'https://www.dropbox.com/s/7ek4lwwf2b7f749/ldblk_1kg_eas.tar.gz?dl=0',
#         }[wc['popul']]
#     wildcard_constraints:
#         popul='[A-Z]+'
#     log:
#         "logs/download_prscs_ld_reference_{popul}.log"
#     shell:
#         "("
#         "wget -O {output} {params[download_link]}"
#         ") &> {log}"


# rule unpack_prscs_ld_reference_1:
#     # unpack the reference for 1000 Genomes
#     # note: the path has to contain "1kg" or "ukbb", otherwise PRScs crashes
#     input:
#         rules.download_prscs_ld_reference.output
#     output:
#         snpinfo='resources/LD_matrix/prscs/1000G/precomputed/{popul}/ldblk_1kg/snpinfo_1kg_hm3',
#         chr=expand('resources/LD_matrix/prscs/1000G/precomputed/{popul}/ldblk_1kg/ldblk_1kg_chr{chr}.hdf5', chr=range(1,23), allow_missing=True)
#     log:
#         "logs/unpack_prscs_ld_reference_{popul}.log"
#     shell:
#         "("
#         "cd resources/LD_matrix/prscs/1000G/precomputed/{wildcards[popul]} && tar -xvf \"$(basename {input})\" && "
#         "if [ -d ldblk_ukbb ]; then rm -r ldblk_ukbb; fi; "
#         "mv ldblk_ukbb_*/ ./ldblk_ukbb && rm \"$(basename {input})\" "
#         ") &> {log}"


# rule run_prscs_precompld_1kg:
#     # 4.4.2 Run PRScs with precomputed 1000G LD reference panel
#     # runs PRScs for the 1000G populations
#     # TODO: replace hardcoded parameters with configurable ones
#     input:
#         super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
#         qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
#         ld_reference=lambda wc: expand(rules.unpack_prscs_ld_reference_1.output, popul=studies.ancestry[studies.study_id == wc.study])
#     params:
#         study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc.study].iloc[0],
#         ld_reference_dir = lambda wc, input: '/'.join(input['ld_reference'][0].split('/')[:-1])
#     output:
#         adj_ss=expand('{geno1kg}/Score_files_for_polygenic/PRScs_precompld/{{study}}/1KGPhase3.w_hm3.{{study}}_pst_eff_a1_b0.5_phi{phi}_chr8.txt', geno1kg=config['Geno_1KG_dir'], phi=['1e+00','1e-02','1e-04','1e-06','auto'], chr=range(1,23)),
#         scale=expand('{geno1kg}/Score_files_for_polygenic/PRScs_precompld/{{study}}/1KGPhase3.w_hm3.{{study}}.{superpop}.scale', geno1kg=config['Geno_1KG_dir'], superpop=config['1kg_superpop']),
#         log='{}/Score_files_for_polygenic/PRScs_precompld/{{study}}/1KGPhase3.w_hm3.{{study}}.log'.format(config['Geno_1KG_dir']),
#         # touch('{}/Score_files_for_polygenic/PRScs_precompld/{{study}}/prscs.ok'.format(config['Geno_1KG_dir']))
#     log:
#         "logs/prs_PRScs/prscs_precompld_1kg_{study}.log"
#     conda:
#         "../envs/py2.yaml"
#     threads:
#         16
#     shell:
#         "( "
#         "{config[Rscript]} {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_PRScs/polygenic_score_file_creator_PRScs.R "
#         "--ref_plink_chr resources/1kg/1KGPhase3.w_hm3.chr "
#         "--ref_keep resources/1kg/keep_files/{params[study_ancestry]}_samples.keep "
#         "--sumstats {input[qc_stats]} "
#         "--plink {config[plink1_9]} "
#         "--memory 5000 "
#         "--output resources/1kg/Score_files_for_polygenic/PRScs_precompld/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
#         "--ref_pop_scale {input[super_pop_keep]} "
#         "--PRScs_path {config[PRScs_dir]} "
#         "--PRScs_ref_path {params[ld_reference_dir]} "
#         "--n_cores {threads} "
#         "--phi_param 1e-6,1e-4,1e-2,1,auto "
#         ") &> {log} "
