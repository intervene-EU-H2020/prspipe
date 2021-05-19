

# https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html

rule harmonize_ukbb:
    # 1.1 UK Biobank
    input:
        rules.extract_hm3.output
    output:
        log='{}/Genotype/Harmonised/UKBB.harmonisation.chr{{chr}}.log'.format(config['UKBB_output']),
        bim='{}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr{{chr}}.bim'.format(config['UKBB_output']),
        bed='{}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr{{chr}}.bed'.format(config['UKBB_output']),
        fam='{}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr{{chr}}.fam'.format(config['UKBB_output'])
    log:
        "logs/harmonize_ukbb_chr{chr}.log"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Harmonisation_of_UKBB/Harmonisation_of_UKBB.R "
        "--chr {wildcards[chr]} "
        "--input_dir {config[UKBB_original]} "
        "--target_fam {config[UKBB_fam]} "
        "--output_dir {config[UKBB_output]}/Genotype/Harmonised "
        "--reference_dir {config[Geno_1KG_dir]} "
        "--qctool2 {config[qctool2]} "
        "--plink {config[plink1_9]} "
        ") &> {log} "


rule ancestry_scoring_ukbb:
    # 2.1 UK Biobank
    # recommended 15G of RAM
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        harmonised_geno=expand(rules.harmonize_ukbb.output, chr=range(1,23))
    output:
        # touch("{}/Projected_PCs/Ancestry_idenitfier/test.ok".format(config['UKBB_output'])),
        eigenvec_sp=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.{superpop}.eigenvec', ukbb_output=config['UKBB_output'], superpop=config['1kg_superpop']),
        keep_sp=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.{superpop}.keep', ukbb_output=config['UKBB_output'], superpop=config['1kg_superpop']),
        scale_sp=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.{superpop}.scale', ukbb_output=config['UKBB_output'], superpop=config['1kg_superpop']),
        eigenvec=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.eigenvec', ukbb_output=config['UKBB_output']),
        eigenvec_var=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.eigenvec.var', ukbb_output=config['UKBB_output']),
        log=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.log', ukbb_output=config['UKBB_output']),
        model_pred_keep=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.model_pred.{superpop}.keep', ukbb_output=config['UKBB_output'], superpop=config['1kg_superpop']),
        model_pred=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.model_pred', ukbb_output=config['UKBB_output']),
        pop_model=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.pop_model.rds', ukbb_output=config['UKBB_output']),
        model_details=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.pop_model_prediction_details.txt', ukbb_output=config['UKBB_output']),
        scale=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry.scale', ukbb_output=config['UKBB_output'])
    log:
        "logs/ancestry_scoring_ukbb.log"
    shell:
        "( "
        "Rscript {config[GenoPred_dir]}/Scripts/Ancestry_identifier/Ancestry_identifier.R "
        "--target_plink_chr {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--n_pcs 100 "
        "--plink {config[plink1_9]} "
        "--plink2 {config[plink2]} "
        "--output {config[UKBB_output]}/Projected_PCs/Ancestry_idenitfier/UKBB.w_hm3.AllAncestry "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--pop_data {config[Geno_1KG_dir]}/integrated_call_samples_v3.20130502.ALL.panel_small "
        ") &> {log} "
        
        
rule ancestry_scoring_ukbb_stringent:
    # get ancestry keep files using a more stringent cutoff of 0.995
    input:
        rules.ancestry_scoring_ukbb.output['model_pred']
    output:
        keep=expand('{ukbb_output}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.{superpop}.keep', ukbb_output=config['UKBB_output'], superpop=config['1kg_superpop'])
    params:
        ancestries=" ".join(sorted(config['1kg_superpop']))
    log:
        'logs/ancestry_scoring_ukbb_stringent.log'
    shell:
        "("
        "bash workflow/scripts/bash/extract_ancestry.sh {input} \"{params[ancestries]}\" {config[UKBB_output]}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry "
        ") &> {log}"
        
        
rule calculate_maf_ancestry_ukbb:
    # calculate allele frequencies for different superpopulations
    input:
        rules.ancestry_scoring_ukbb_stringent.output,
        rules.harmonize_ukbb.output
    output:
        afreq_gz='{}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.{{superpop}}.chr{{chr}}.afreq.gz'.format(config['UKBB_output'])
    log:
        'logs/calculate_maf_ancestry_ukbb_{superpop}_{chr}.log'
    shell:
        '('
        '{config[plink2]} '
        '--bfile {config[UKBB_output]}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr{wildcards[chr]} '
        '--freq '
        '--keep {config[UKBB_output]}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.{wildcards[superpop]}.keep '
        '--out {config[UKBB_output]}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.{wildcards[superpop]}.chr{wildcards[chr]}; '
        'gzip {config[UKBB_output]}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.{wildcards[superpop]}.chr{wildcards[chr]}.afreq; '
        ') &> {log} '
    

rule all_calculate_maf_ancestry_ukbb:
    # runs rule above for all ancestries and chromosomes
    input:
        expand(rules.calculate_maf_ancestry_ukbb.output, chr=range(1,23), superpop=config['1kg_superpop'])
    output:
        touch('{}/Projected_PCs/Ancestry_idenitfier/stringent/UKBB.w_hm3.AllAncestry.all_afreq.ok'.format(config['UKBB_output']))