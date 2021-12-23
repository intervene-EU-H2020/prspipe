# An example collection of Snakemake rules imported in the main Snakefile.

rule download_test_data:
    # rule to download test data from figshare
    shell:
        "wget -O {config[Geno_1KG_dir]}/PRS.tar.gz https://figshare.com/ndownloader/files/30799843?private_link=94afe1b02cf622234566 && "
        "cd {config[Geno_1KG_dir]} && "
        "tar -xvf PRS.tar.gz && "
        "find ./Score_files_for_polygenic -type f -exec touch {{}} +" # touching the extracted files so the snakemake timestamp checks work


rule install_software:
    # in case ./install_software.sh has not been run yet,
    # this rule can be used to do it automatically
    # this should only be used for binaries or software that can not
    # be downloaded through conda, or where using docker/singularity is not necessary
    output:
        "bin/plink",
        "bin/plink2",
        "bin/gcta/gcta64",
        "bin/gctb/gctb"
    shell:
        "bash install_software.sh"


rule download_integrated_call_samples_v3:
    output:
         '{}/integrated_call_samples_v3.20130502.ALL.panel'.format(config['Geno_1KG_dir']),
         '{}/integrated_call_samples_v3.20130502.ALL.panel_small'.format(config['Geno_1KG_dir'])
    shell:
        "("
        "mkdir -p {config[Geno_1KG_dir]} && "
        "cd {config[Geno_1KG_dir]} && "
        "wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel && "
        "cut -f 1-3 integrated_call_samples_v3.20130502.ALL.panel > integrated_call_samples_v3.20130502.ALL.panel_small "
        ")"


rule create_ancestry:
    # 2.3 1000 Genomes populations
    # https://opain.github.io/GenoPred/Pipeline_prep.html
    # Here we download information on which population each individual in the 1000 Genomes reference is from. Individuals are grouped into ‘populations’ which are typically country specific, and ‘super populations’ which include a collection of ancetrally similar countries. Individuals are grouped into these populations if the last few generations of their family are all from one region. We need this population data so we can select individuals in the 1000 Genomes data to match those in our target samples. This is important for providing accurate information on LD structure and minor allele frequencies.
    # TODO: figure out how to handle logging when using "script" directive for R scripts
    input:
        rules.download_integrated_call_samples_v3.output
    output:
        super_pop_keep="{}/super_pop_keep.list".format(config['Geno_1KG_dir']),
        pop_keep="{}/pop_keep.list".format(config['Geno_1KG_dir']),
        super_pop_and_pop="{}/super_pop_and_pop_keep.list".format(config['Geno_1KG_dir']),
        create_ancestry_ok=touch("{}/keep_files/create_ancestry.ok".format(config['Geno_1KG_dir'])),
        keep_files=expand("{}/keep_files/{{popul}}_samples.keep".format(config['Geno_1KG_dir']),popul=config['1kg_superpop'])
    singularity:
        config['singularity']['all']
    resources:
    	threads=1,
        mem_mb=4000,
        misc='--container-image=/dhc/projects/intervene/prspipe_0_0_1.sqsh',
        time="04:00:00"
    script:
        "../scripts/R/setup/create_ancestry.R"


rule download_hapmap3_snplist:
    input:
        'resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping.tsv.gz'
    output:
        "{}/w_hm3.snplist".format(config['HapMap3_snplist_dir'])
    log:
        "logs/download_hapmap3_snplist.log"
    shell:
        "("
        "zcat {input} | cut -f2-4 | "
        ' awk \'BEGIN{{OFS="\t"; print "SNP", "A1", "A2"}}{{if(NR > 1){{print $0}}}}\' > {output} '
        ") &> {log}"
        
        # old version: use the LDSC list
        # "cd {config[HapMap3_snplist_dir]} && "
        # "wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 && "
        # "bunzip2 w_hm3.snplist.bz2"

        
        
rule cleanup_after_setup:
    # removes unnecessary intermediate output files
    log:
        "logs/cleanup.log"
    shell:
        "("
        "rm -f {config[Geno_1KG_dir]}/ALL.chr*.phase3_shapeit2_mvncall_integrated_v*.*.genotypes.vcf.gz ; "
        "rm -f {config[Geno_1KG_dir]}/merge_list.txt ; "
        "rm -f {config[Geno_1KG_dir]}/1KGPhase3.chr*.bed ; "
        "rm -f {config[Geno_1KG_dir]}/1KGPhase3.chr*.bim ; "
        "rm -f {config[Geno_1KG_dir]}/1KGPhase3.chr*.fam ; "
        "rm -f {config[Geno_1KG_dir]}/1KGPhase3.chr*.nosex ; "
        "rm -f {config[Geno_1KG_dir]}/1KGPhase3.chr*.log ; "
        "rm -f {config[Geno_1KG_dir]}/1KGPhase3.chr*.extract ; "
        "rm -f {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr*.nosex ; "
        "rm -f {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr*.log ; "
        "rm -f {config[Geno_1KG_dir]}/freq_files/*/*log ; "
        "rm -f {config[Geno_1KG_dir]}/freq_files/*/*noseq ; "
        ") &> {log}"
