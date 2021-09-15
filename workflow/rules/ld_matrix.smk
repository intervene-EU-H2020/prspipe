# Rules for downloading LD references or calculating them from scratch

rule create_dbslmm_sblup_ld_ref:
    # Compute a custom LD reference panel 
    # This is computed from the .bed/.bim/.fam files with the prefix given by the bfile input (here, 1000G data)
    output:
        ld_ref="{}/sblup_dbslmm/1000G/fromscratch/{{ancestry}}/{{chromosome}}.l2.ldscore.gz".format(config['LD_ref_dir'])
    log:
        "logs/ld_matrix/dbslmm.{ancestry}.{chromosome}.log"
    conda:
        "../envs/ldsc.yaml"
    shell:
        "python {config[LDSC_dir]}/ldsc.py "
        "--bfile {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr{wildcards[chromosome]} "
        "--keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--l2 "
        "--ld-wind-cm 1 "
        "--yes-really "
        "--out {config[LD_ref_dir]}/sblup_dbslmm/1000G/fromscratch/{wildcards[ancestry]}/{wildcards[chromosome]}"


rule download_dbslmm_sblup_ld_ref:
    # Download a (precomputed) LD reference panel for the given ancestry based on 1000G (note that it currently only supports EUR)
    output:
        ld_ref=expand("{ldrefdir}/sblup_dbslmm/1000G/precomputed/{{ancestry}}/{chromosome}.l2.ldscore.gz", ldrefdir=config['LD_ref_dir'], chromosome=range(1,23))
    shell:
        "mkdir -p {config[LD_ref_dir]}/sblup_dbslmm/1000G/precomputed/{wildcards[ancestry]} ;"
        "cd {config[LD_ref_dir]}/sblup_dbslmm/1000G/precomputed/EUR ;"
        "wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 ;"
        "tar -xf eur_w_ld_chr.tar.bz2 && rm eur_w_ld_chr.tar.bz2 ;"
        "mv eur_w_ld_chr/* . && rm -r eur_w_ld_chr"

