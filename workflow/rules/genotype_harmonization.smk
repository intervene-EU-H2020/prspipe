
#########################################
# harmonization of target genotype data #
#########################################

def get_input_file_from_type(prefix, chr, type):
    # get an input file to check for based on prefix, chromosome, and input file type
    if type == 'samp_imp_plink1':
        return(prefix + str(chr) + '.bim')
    elif type == 'samp_imp_vcf':
        return(prefix + str(chr) + '.vcf')
    else:
        if 'samp_imp_bgen' in type:
            return(prefix + str(chr) + '.bgen')
        else:
            raise NotImplementedError('Error: type {} not implemented (check your target_list: {})'.format(type, config['target_list']))

def type_to_arg_geno_to_plink(type):
    # return the correct arguments for geno_to_plink.R
    try:
        return {
            'samp_imp_bgen_ref_first':('samp_imp_bgen', 'ref-first'),
            'samp_imp_bgen_ref_last':('samp_imp_bgen', 'ref-last'),
            'samp_imp_bgen':('samp_imp_bgen','ref-unknown')}[type]
    except KeyError:
        return (type, 'ref-last')


rule harmonize_target_genotypes:
    # target sample harmonization using Ollie's script
    input:
        bim=rules.extract_hm3.output['bim'],
        fam=rules.extract_hm3.output['fam'],
        bed=rules.extract_hm3.output['bed'],
        geno = lambda wc: get_input_file_from_type(target_list.loc[wc['bbid'], 'path'] if wc['bbid'] in target_list.index else '',
                                                   wc['chr'],
                                                   target_list.loc[wc['bbid'], 'type'] if wc['bbid'] in target_list.index else 'samp_imp_plink1'),
        liftover=rules.download_liftover.output['liftover'],
        chain=rules.download_liftover.output['hg19_to_hg38_chain']
    output:
        bim='custom_input/{bbid}/genotypes/chr{chr}.bim',
        bed='custom_input/{bbid}/genotypes/chr{chr}.bed',
        fam='custom_input/{bbid}/genotypes/chr{chr}.fam'
    params:
        in_prefix_target = lambda wc, input: ''.join(input['geno'].split('.')[:-1]),
        in_prefix_ref = lambda wc, input: input['bim'][:-4],
        out_prefix=lambda wc, output: output['bim'][:-4],
        arg = lambda wc: type_to_arg_geno_to_plink(target_list.loc[wc['bbid'], 'type'])
    threads:
        4
    resources:
        mem_mb=10000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home",
        time='12:00:00'
    log:
        'logs/harmonize_target_genotypes/{bbid}/{chr}.log'
    shell:
        '('
        '{config[Rscript]} workflow/scripts/GenoPred/Scripts/geno_to_plink/geno_to_plink.R '
        '--target {params[in_prefix_target]} '
        '--ref {params[in_prefix_ref]} '
        '--format {params[arg][0]} '
        '--plink {config[plink1_9]} '
        '--plink2 {config[plink2]} '
        '--liftover {input[liftover]} '
        '--liftover_track {input[chain]} '
        '--out {params[out_prefix]} '
        '--bgen_ref {params[arg][1]} '
        '--threads {threads} '
        '--mem_mb {resources[mem_mb]} '
        ') &> {log} '
        
        
rule all_harmonize_target_genotypes:
    input:
        expand(rules.harmonize_target_genotypes.output, bbid=target_list.index.tolist(), chr=range(1,23))
