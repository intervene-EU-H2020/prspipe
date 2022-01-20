
# Rules for downloading the shared LD references for DBSLMM / SBLUP 

rule create_dbslmm_sblup_ld_ref:
    # Compute a custom LD reference panel 
    # This is computed from the .bed/.bim/.fam files with the prefix given by the bfile input (here, 1000G data)
    output:
        ld_ref="resources/LD_matrix/sblup_dbslmm/1000G/fromscratch/{ancestry}/{chr}.l2.ldscore.gz"
    log:
        "logs/create_dbslmm_sblup_ld_ref/{ancestry}.{chr}.log"
    conda:
        "../envs/ldsc.yaml"
    shell:
        "python {config[LDSC_dir]}/ldsc.py "
        "--bfile resources/1kg/1KGPhase3.w_hm3.chr{wildcards[chr]} "
        "--keep resources/1kg/keep_files/{wildcards[ancestry]}_samples.keep "
        "--l2 "
        "--ld-wind-cm 1 "
        "--yes-really "
        "--out resources/1kg/sblup_dbslmm/1000G/fromscratch/{wildcards[ancestry]}/{wildcards[chr]}"


rule download_dbslmm_sblup_ld_ref:
    # Download a (precomputed) LD reference panel for the given ancestry based on 1000G (note that it currently only supports EUR)
    output:
        ld_ref=expand("resources/LD_matrix/sblup_dbslmm/1000G/precomputed/{{ancestry}}/{chr}.l2.ldscore.gz", chr=range(1,23))
    shell:
        "mkdir -p resources/LD_matrix/sblup_dbslmm/1000G/precomputed/{wildcards[ancestry]} ;"
        "cd resources/LD_matrix/sblup_dbslmm/1000G/precomputed/EUR ;"
        "wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 ;"
        "tar -xf eur_w_ld_chr.tar.bz2 && rm eur_w_ld_chr.tar.bz2 ;"
        "mv eur_w_ld_chr/* . && rm -r eur_w_ld_chr"

