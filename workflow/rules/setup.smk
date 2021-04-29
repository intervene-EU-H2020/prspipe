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
    # TODO: figure out how to handle logging when using "script" directive for R scripts
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
    log:
        "logs/get_plink_files_chr/{chr}.log"
    shell:
        # v5b version does not contain rs ids!
        # "(cd {config[Geno_1KG_dir]} && wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{wildcards[chr]}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz); "
        # "{config[plink1_9]} --vcf {config[Geno_1KG_dir]}/ALL.chr{wildcards[chr]}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz "
        "(cd {config[Geno_1KG_dir]} && wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140708_previous_phase3/v5_vcfs/ALL.chr{wildcards[chr]}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz) 2> {log}; "
        "({config[plink1_9]} --vcf {config[Geno_1KG_dir]}/ALL.chr{wildcards[chr]}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz "
        "--make-bed "
        "--out {config[Geno_1KG_dir]}/1KGPhase3.chr{wildcards[chr]};) &> {log} "
        "# rm {config[Geno_1KG_dir]}/ALL.chr{wildcards[chr]}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz "


rule get_plink_files_chr_all:
    # run rule above for all chromosomes
    input:
        expand(rules.get_plink_files_chr.output, chr=range(1,23))


rule download_ensembl_variation_vcf:
    # currently not used.
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
        # ensembl_vcf=rules.download_ensembl_variation_vcf.output
    output:
        extract=expand("{}/1KGPhase3.chr{{chr}}.extract".format(config['Geno_1KG_dir']), chr=range(1,23))
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
    log:
        "logs/extract_hm3/extract_hm3.log"
    shell:
        "("
        "for chr in $(seq 1 22); do "
        "{config[plink1_9]} --bfile {config[Geno_1KG_dir]}/1KGPhase3.chr${{chr}} "
        "--make-bed "
        "--extract {config[Geno_1KG_dir]}/1KGPhase3.chr${{chr}}.extract "
        "--out {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr${{chr}}; " 
        "done; "
        "ls {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr*.bed | sed -e 's/\.bed//g' > {config[Geno_1KG_dir]}/merge_list.txt ;"
        "{config[plink1_9]} --merge-list {config[Geno_1KG_dir]}/merge_list.txt "
        "--make-bed "
        "--out {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        ") &> {log} "

rule allele_freq_pop:
    # 2.5 1000 Genomes allele frequency files
    # Here we create files containing ancestry specific minor allele frequency estimates for the HapMap3 SNPs based on the 1000 Genomes Phase 3 data. This information is mainly used for mean-imputation of missing SNPs during genotype-based scoring. This avoids target sample specific minor allele frequencies being used for mean imputation which may not be available, and will vary between target samples.
    input:
        rules.create_ancestry.output,
        rules.extract_hm3.output
    params:
        keep_file=lambda wc: "{}/keep_files/{}_samples.keep".format(config['Geno_1KG_dir'], wc['popul'])
    output:
        pop_list=expand("{}/freq_files/{{popul}}/1KGPhase3.w_hm3.{{popul}}.chr{{chr}}.frq".format(config['Geno_1KG_dir']), chr=range(1,23), allow_missing=True)
    wildcard_constraints:
        pop='[A-Z]+'
    log:
        "logs/allele_freq_pop/{popul}.log"
    shell:
        "("
        "for chrom in $(seq 1 22); do "
        "{config[plink1_9]} --bfile {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr${{chrom}} "
        "--freq "
        "--out {config[Geno_1KG_dir]}/freq_files/{wildcards[popul]}/1KGPhase3.w_hm3.{wildcards[popul]}.chr${{chrom}}; "
        "done "
        ") &> {log} "

rule run_allele_freq_pop:
    # runs 2.5 for populations
    input:
        rules.create_ancestry.output,
        expand(rules.allele_freq_pop.output, popul=config['1kg_pop'])
    output:
        touch("{}/freq_files/all_pop.ok".format(config['Geno_1KG_dir']))


rule run_allele_freq_superpop:
    # runs 2.5 for super populations
    input:
        rules.create_ancestry.output,
        expand(rules.allele_freq_pop.output, popul=config['1kg_superpop'])
    output:
        touch("{}/freq_files/all_superpop.ok".format(config['Geno_1KG_dir']))


rule run_allele_freq_allancestry:
    # runs 2.5 for all ancestries combined
    input:
        rules.extract_hm3.output
    output:
        expand("{}/freq_files/AllAncestry/1KGPhase3.w_hm3.AllAncestry.chr{{chr}}.frq".format(config['Geno_1KG_dir']), chr=range(1, 23))
    log:
        "logs/run_allele_freq_allancestry/run_allele_freq_allancestry.log"
    shell:
        "("
        "for chr in $(seq 1 22); do "
        "{config[plink1_9]} --bfile {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr${{chr}} "
        "--freq "
        "--out {config[Geno_1KG_dir]}/freq_files/AllAncestry/1KGPhase3.w_hm3.AllAncestry.chr${{chr}}; "
        "done "
        ") &> {log} "


rule run_allele_freq_all:
    # this rule executes all the rules above
    input:
        rules.run_allele_freq_pop.output,
        rules.run_allele_freq_superpop.output,
        rules.run_allele_freq_allancestry.output
    input:
        touch("{}/freq_files/all.ok".format(config['Geno_1KG_dir']))


rule ancestry_scoring:
    # 3 Ancestry scoring
    # https://opain.github.io/GenoPred/Pipeline_prep.html
    # In this section we perform principal components analysis of genotypic data in the 1000 Genomes reference to identify the main axes of population structure, which typically correspond to ancestral differences. We calculate these principal components (PCs) of ancestry across the full reference, and within super populations to detect broad and ancestry-specific axes of variance. After PCA, we idenitfy which variants are associated with the PCs (SNP-weights) to calculate ancestry scores on the same axes of variance in future target samples. This can be used to infer the ancestry of an individual which is an important factor to consider when performing genotype-baed prediction.
    #TODO: add depedency handling: install.packages(c('data.table','caret','pROC','verification','optparse'))
    # The script will always create .eigenvec and .eigenvec.var files. These contain the PCs score for individuals in the reference dataset, and the SNP-weights for the PCs respectively.
    input:
        rules.create_ancestry.output,
        rules.extract_hm3.output
    output:
        # pop_enet_model="{}/Score_files_for_ancestry/{{popul}}/1KGPhase3.w_hm3.{{popul}}.pop_enet_model.rds".format(config['Geno_1KG_dir']),
        eigenvec="{}/Score_files_for_ancestry/{{popul}}/1KGPhase3.w_hm3.{{popul}}.eigenvec".format(config['Geno_1KG_dir']),
        pc_scale="{}/Score_files_for_ancestry/{{popul}}/1KGPhase3.w_hm3.{{popul}}.scale".format(config['Geno_1KG_dir']),
        eigenvec_var="{}/Score_files_for_ancestry/{{popul}}/1KGPhase3.w_hm3.{{popul}}.eigenvec.var".format(config['Geno_1KG_dir'])
    log:
        # the script is configured to write log files here
        "{}/Score_files_for_ancestry/{{popul}}/1KGPhase3.w_hm3.{{popul}}.log".format(config['Geno_1KG_dir'])
    wildcard_constraints:
        popul='[A-Z]+'
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/ancestry_score_file_creator/ancestry_score_file_creator.R "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[popul]}_samples.keep "
        "--plink {config[plink1_9]} "
        "--plink2 {config[plink2]} "
        "--n_pcs 100 "
        "--output {config[Geno_1KG_dir]}/Score_files_for_ancestry/{wildcards[popul]}/1KGPhase3.w_hm3.{wildcards[popul]}"
        ") &> {log}"


rule ancestry_scoring_allancestry:
    # 3 Ancestry scoring
    # https://opain.github.io/GenoPred/Pipeline_prep.html
    # If --ref_pop_scale is specified, the script will also creates files stating the mean and standard deviation of the PCs for each group. Furthermore, it will derive an elastic net model predicting each group, and report the accuracy of the derived models.
    input:
        rules.create_ancestry.output,
        rules.extract_hm3.output
    output:
        pop_enet_model="{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.pop_enet_model.rds".format(config['Geno_1KG_dir']),
        eigenvec="{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.eigenvec".format(config['Geno_1KG_dir']),
        pc_scale="{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.scale".format(config['Geno_1KG_dir']),
        eigenvec_var="{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.eigenvec.var".format(config['Geno_1KG_dir']),
    log:
        # the script is configured to write log files here
        "{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.log".format(config['Geno_1KG_dir'])
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/ancestry_score_file_creator/ancestry_score_file_creator.R "
        "--ref_plink_chr {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr "
        "--plink {config[plink1_9]} "
        "--plink2 {config[plink2]} "
        "--n_pcs 100 "
        "--ref_pop_scale {config[Geno_1KG_dir]}/super_pop_keep.list "
        "--output {config[Geno_1KG_dir]}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry"
        ") &> {log}"


rule all_setup:
    # this triggers all the steps 1-3 (all "setup" steps)
    # the Polygenic scoring steps (>=4) that follow are implemented in a different module
    input:
        '{}/Score_files_for_ancestry/EUR/1KGPhase3.w_hm3.EUR.eigenvec.var'.format(config['Geno_1KG_dir']),
        '{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.eigenvec.var'.format(config['Geno_1KG_dir']),
        rules.run_allele_freq_all.output
