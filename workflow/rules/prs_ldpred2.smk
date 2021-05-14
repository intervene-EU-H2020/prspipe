
# TODO: add dependencies bigsnpr, bigreadr, runonce

rule generate_ld_referece_ldpred2:
    # 4.8.1 Create LD reference for LDPred2, part 1
    # uses a lot of memory > 16G
    # Mostly single core seems to be used, adding cores probably won't speed it up much
    # TODO: move parameters to config (?)
    input:
        keep=rules.create_ancestry.output['super_pop_keep'],
        genetic_map = rules.download_genetic_map.output,
        bed_bim_fam = rules.extract_hm3.output
    output:
        rds=expand('{geno1kdir}/LD_matrix/LDPred2/{popul}/LD_chr{chr}.rds', geno1kdir=config['Geno_1KG_dir'], chr=range(1,23), allow_missing=True),
        map=expand('{geno1kdir}/LD_matrix/LDPred2/{popul}/map_tmp.rds', geno1kdir=config['Geno_1KG_dir'], allow_missing=True),
        sd=expand('{geno1kdir}/LD_matrix/LDPred2/{popul}/sd_tmp.rds', geno1kdir=config['Geno_1KG_dir'], allow_missing=True)
    params:
        cores=6
    log:
        "logs/generate_ld_referece_ldpred2_{popul}.log"
    shell:
        "("
        "Rscript workflow/scripts/R/prs_ldpred2/generate_ld_reference.R {wildcards[popul]} {params[cores]} "
        ") &> {log} "
        
        
rule download_ld_reference_ldpred2:
    # download EUR (UKBB) LD reference
    output:
        touch('{}/LD_matrix/LDPred2_precomputed/EUR/dl_ld.txt'.format(config['Geno_1KG_dir']))
    shell:
        "cd \"$(dirname {output})\"; "
        "wget -O ldpred2_reference.zip https://ndownloader.figshare.com/articles/13034123/versions/3 ; "
        "unzip ldpred2_reference.zip; "
        

#rule generate_sd_file_ldpred2:
#    # 4.8.1 Create LD reference for LDPred2, part 2
#    # RM: not needed. 
#    input:
#        map=rules.generate_ld_referece_ldpred2.output.map
#    output:
#        sd='{}/LD_matrix/LDPred2/{{popul}}/sd.rds'.format(config['Geno_1KG_dir']),
#        map='{}/LD_matrix/LDPred2/{{popul}}/map.rds'.format(config['Geno_1KG_dir'])
#    log:
#        "logs/generate_sd_file_ldpred2_{popul}.log"
#    shell:
#        "("
#        "Rscript workflow/scripts/R/prs_ldpred2/generate_sd_file.R {input[map]} " 
#        ") &> {log}"


