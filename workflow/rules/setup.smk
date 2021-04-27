# An example collection of Snakemake rules imported in the main Snakefile.

rule install_software:
    # in case ./install_software.sh has not been run yet,
    # this rule can be used to do it automatically
    # this should only be used for binaries or software that can not
    # be downloaded through conda
    output:
        "bin/plink",
        "bin/plink2",
        "bin/gcta/gcta64",
        "bin/gctb/gctb"
    shell:
        "bash install_software.sh"


rule create_ancestry:
    # 2.3 1000 Genomes populations
    # https://opain.github.io/GenoPred/Pipeline_prep.html
    # Here we download information on which population each individual in the 1000 Genomes reference is from. Individuals are grouped into ‘populations’ which are typically country specific, and ‘super populations’ which include a collection of ancetrally similar countries. Individuals are grouped into these populations if the last few generations of their family are all from one region. We need this population data so we can select individuals in the 1000 Genomes data to match those in our target samples. This is important for providing accurate information on LD structure and minor allele frequencies.
    output:
        super_pop_keep="{}/super_pop_keep.list".format(config['Geno_1KG_dir']),
        pop_keep="{}/pop_keep.list".format(config['Geno_1KG_dir']),
        super_pop_and_pop="{}/super_pop_and_pop_keep.list".format(config['Geno_1KG_dir']),
        create_ancestry_ok=touch("{}/keep_files/create_ancestry.ok".format(config['Geno_1KG_dir']))
    script:
        "../scripts/R/setup/create_ancestry.R"


rule get_plink_files_chr:
    # 2.4 1000 Genomes PLINK files (step 1)
    # https://opain.github.io/GenoPred/Pipeline_prep.html
    # Download vcf file, convert to plink
    output:
         bed="{}/1KGPhase3.chr{{chr}}.bed".format(config['Geno_1KG_dir']),
         fam="{}/1KGPhase3.chr{{chr}}.fam".format(config['Geno_1KG_dir']),
         bim="{}/1KGPhase3.chr{{chr}}.bim".format(config['Geno_1KG_dir'])
    shell:
        # v5b version does not contain rs ids!
        # "(cd {config[Geno_1KG_dir]} && wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{wildcards[chr]}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz); "
        # "{config[plink1_9]} --vcf {config[Geno_1KG_dir]}/ALL.chr{wildcards[chr]}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz "
        "(cd {config[Geno_1KG_dir]} && wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140708_previous_phase3/v5_vcfs/ALL.chr{wildcards[chr]}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz); "
        "{config[plink1_9]} --vcf {config[Geno_1KG_dir]}/ALL.chr{wildcards[chr]}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz "
        "--make-bed "
        "--out {config[Geno_1KG_dir]}/1KGPhase3.chr{wildcards[chr]}; "
        "# rm {config[Geno_1KG_dir]}/ALL.chr{wildcards[chr]}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz "


rule get_plink_files_chr_all:
    # run rule above for all chromosomes
    input:
        expand(rules.get_plink_files_chr.output, chr=range(1,23))


rule download_ensembl_variation_vcf:
    # not needed when downloading the archived v5!
    # download the 1KG vcf from Ensembl, in order to get the rs ids
    #TODO add requirement bcftools 
    #TODO remove temporary output files
    output:
        expand("{}/1000GENOMES-phase_3.chr{{chr}}.vcf.gz".format(config['Geno_1KG_dir']), chr=range(1,23))
    shell:
        "cd {config[Geno_1KG_dir]}; "
        "if [ ! -f 1000GENOMES-phase_3.vcf.gz ]; then "
        "wget http://ftp.ensembl.org/pub/grch37/release-103/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz && "
        "wget http://ftp.ensembl.org/pub/grch37/release-103/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz.csi; "
        "fi;"
        "for chr in $(seq 1 22); do "
        "bcftools view -r ${{chr}} -O z 1000GENOMES-phase_3.vcf.gz > 1000GENOMES-phase_3.chr${{chr}}.vcf.gz ; "
        "done"
        

rule snp_to_iupac:
    # 2.4 1000 Genomes PLINK files (step 2)
    # https://opain.github.io/GenoPred/Pipeline_prep.html
    # Use R to identify a list of SNPs that matches with the hapmap3 snplist
    input:
        bed_bim_fam=rules.get_plink_files_chr_all.input,
        hapmap3_snplist='{}/w_hm3.snplist'.format(config['HapMap3_snplist_dir']),
        ensembl_vcf=rules.download_ensembl_variation_vcf.output
    output:
        extract=expand("{}/1KGPhase3.chr{{chr}}.extract".format(config['Geno_1KG_dir']), chr=range(19,21))
    script:
        "../scripts/R/setup/snp_to_iupac.R"


rule extract_hm3:
    # 2.4 1000 Genomes PLINK files (step 3)
    # https://opain.github.io/GenoPred/Pipeline_prep.html
    # extract HapMap 3 SNPs and merge into single .bed
    # TODO: intermediate outputs from rule "get_plink_files_chr_all" can be removed after.
    input:
        rules.snp_to_iupac.output
    output:
        bed='{}/1KGPhase3.w_hm3.GW.bed'.format(config['Geno_1KG_dir']),
        bim='{}/1KGPhase3.w_hm3.GW.bim'.format(config['Geno_1KG_dir']),
        fam='{}/1KGPhase3.w_hm3.GW.fam'.format(config['Geno_1KG_dir'])
    shell:
        "for chr in $(seq 19 20); do "
        "{config[plink1_9]} --bfile {config[Geno_1KG_dir]}/1KGPhase3.chr${{chr}} "
        "--make-bed "
        "--extract {config[Geno_1KG_dir]}/1KGPhase3.chr${{chr}}.extract "
        "--out {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr${{chr}}; " 
        "done; "
        "ls {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr*.bed | sed -e 's/\.bed//g' > {config[Geno_1KG_dir]}/merge_list.txt ;"
        "{config[plink1_9]} --merge-list {config[Geno_1KG_dir]}/merge_list.txt "
        "--make-bed "
        "--out {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "

