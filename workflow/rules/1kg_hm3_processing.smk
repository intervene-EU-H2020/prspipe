
##################################################
# Download and process the 1000Genomes reference #
##################################################

if 'wget_use_https' in config:
    if config['wget_use_https'] == True :
        ftp_prefix='https'
    elif config['wget_use_https'] == False :
        ftp_prefix='ftp'
    else:
        raise ValueError('wget_use_https has to be True or False, if specified. Got "{}"'.format(config['wget_use_https']))
else:
    ftp_prefix='ftp'


rule download_1kg:
    # generate 1000 Genomes PLINK files
    # Download vcf file, convert to plink
    output:
        bim="resources/1kg/1KGPhase3.chr{chr}.bim",
        bed="resources/1kg/1KGPhase3.chr{chr}.bed",
        fam="resources/1kg/1KGPhase3.chr{chr}.fam"
    params:
        # v5b - no rsIDs!: ftp_path=lambda wc: "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz".format(wc['chr']),
        # v5 - has rsIDs:
        ftp_path=lambda wc: "{}://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140708_previous_phase3/v5_vcfs/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz".format(ftp_prefix, wc['chr']),
        out = lambda wc, output: output['bed'][:-4]
    log:
        "logs/download_1kg/{chr}.log"
    singularity:
        config['singularity']['all']
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    shell:
        "("
        "export TMPDIR=\"$(readlink -f ./temp)\"; "
        "if [ -f resources/1kg/chr{wildcards[chr]}.vcf.gz ]; then rm resources/1kg/chr{wildcards[chr]}.vcf.gz; fi; "
        "if [ -f resources/1kg/chr{wildcards[chr]}.vcf ]; then rm resources/1kg/chr{wildcards[chr]}.vcf; fi; "
        "wget {params[ftp_path]} -nv -O resources/1kg/chr{wildcards[chr]}.vcf.gz && "
        "{config[plink1_9]} --vcf resources/1kg/chr{wildcards[chr]}.vcf.gz "
        "--make-bed "
        "--out {params[out]} && "
        "rm resources/1kg/chr{wildcards[chr]}.vcf.gz && "
        "rm {params[out]}.log "
        ") &> {log} "


rule download_integrated_call_samples_v3:
    output:
         'resources/1kg/integrated_call_samples_v3.20130502.ALL.panel',
         'resources/1kg/integrated_call_samples_v3.20130502.ALL.panel_small'
    params:
        url = '{}://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel'.format(ftp_prefix)
    shell:
        "("
        "mkdir -p resources/1kg && "
        "cd resources/1kg && "
        "wget {params[url]} && "
        "cut -f 1-3 integrated_call_samples_v3.20130502.ALL.panel > integrated_call_samples_v3.20130502.ALL.panel_small "
        ")"

rule create_ancestry:
    # 2.3 1000 Genomes populations
    # https://opain.github.io/GenoPred/Pipeline_prep.html
    input:
        rules.download_integrated_call_samples_v3.output
    output:
        super_pop_keep="resources/1kg/super_pop_keep.list",
        pop_keep="resources/1kg/pop_keep.list",
        super_pop_and_pop="resources/1kg/super_pop_and_pop_keep.list",
        create_ancestry_ok=touch("resources/1kg/keep_files/create_ancestry.ok"),
        keep_files=expand("resources/1kg/keep_files/{ancestry}_samples.keep", ancestry=config['1kg_superpop'])
    singularity:
        config['singularity']['all']
    resources:
        threads=1,
        time="04:00:00",
        mem_mb=4000,
        misc='--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home '
    script:
        "../scripts/R/setup/create_ancestry.R"




rule download_liftover:
    # Note: won't be necessary if run with singularity
    output:
        liftover="liftover/liftOver",
        hg18_to_hg19_chain="liftover/hg18ToHg19.over.chain",
        hg18_to_hg38_chain="liftover/hg18ToHg38.over.chain",
        hg19_to_hg38_chain="liftover/hg19ToHg38.over.chain"
    log:
        "logs/download_liftover.log"
    shell:
        "("
        "mkdir -p liftover && "
        "cd liftover && "
        "wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && "
        "chmod u+x liftOver && "
        "wget https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz && "
        "gunzip hg18ToHg19.over.chain.gz && "
        "wget https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz && "
        "gunzip hg18ToHg38.over.chain.gz && "
        "wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz && "
        "gunzip hg19ToHg38.over.chain.gz "
        ") &> {log} "


##############################################
# harmonizing the HapMap3 variants - skipped #
##############################################

# Remo: these rules are here for documentation purposes only
#       these rules will be skipped, because we distibute their output with the pipeline.


rule download_hapmap3_r3:
    # release 3
    # skipped.
    output:
        ped="resources/hapmap3/hapmap3_r3_b36_fwd.consensus.qc.poly.ped",
        map="resources/hapmap3/hapmap3_r3_b36_fwd.consensus.qc.poly.map"
    log:
        "logs/download_hapmap3_r3.log"
    shell:
        "("
        "mkdir -p resources/hapmap3 && "
        "cd resources/hapmap3 && "
        "wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.ped.gz && "
        "wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.map.gz && "
        "gunzip *.gz "
        ") &> {log}"
        
        
rule download_hapmap3_r2:
    # release 2
    # skipped.
    output:
        ped="resources/hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.ped",
        map="resources/hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.map"
    log:
        "logs/download_hapmap3_r2.log"
    shell:
        "("
        "mkdir -p resources/hapmap3 && "
        "cd resources/hapmap3 && "
        "wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes//2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2 && "
        "wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes//2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2 && "
        "bunzip2 *.bz2 "
        ") &> {log}"
        
        
rule download_hapmap3_r1:
    # release 1
    # skipped.
    output:
        ped="resources/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.ped",
        map="resources/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map"
    log:
        "logs/download_hapmap3_r1.log"
    shell:
        "("
        "mkdir -p resources/hapmap3 && "
        "cd resources/hapmap3 && "
        "wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2 && "
        "wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.ped.bz2 && "
        "bunzip2 hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2 && "
        "bunzip2 hapmap3_r1_b36_fwd_consensus.qc.poly.recode.ped.bz2 "
        ") &> {log} "


rule make_hapmap3_bim:
    input:
        ped=[rules.download_hapmap3_r1.output['ped'], rules.download_hapmap3_r2.output['ped'], rules.download_hapmap3_r3.output['ped']],
        map=[rules.download_hapmap3_r1.output['map'], rules.download_hapmap3_r2.output['map'], rules.download_hapmap3_r3.output['map']]
    output:
        bim_r1="resources/hapmap3/hapmap3_r1.bim",
        bim_r2="resources/hapmap3/hapmap3_r2.bim",
        bim_r3="resources/hapmap3/hapmap3_r3.bim"
    log:
        "logs/make_hapmap3_bim.log"
    run:
        for i, (ped, map_) in enumerate(zip(input['ped'], input['map'])):
            shell(config['plink1_9'] + " --ped " + ped + " --map " + map_ + ' --out resources/hapmap3/hapmap3_r' + str(i+1), ' --make-just-bim &> ' + str(log) )

        
# this parameter sets the version of the dbSNP R-package used by the rules below.
# as of 4.11.2021, this was the latest release available for GRCh38

dbsnp_v_GRCh38 = '151'

rule install_bioconductor_dbsnp_GRCh38:
    params:
        # version of dbSNP package to use, see https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP151.GRCh38.html
        dbsnp_version=dbsnp_v_GRCh38
    output:
        touch('R/install_dbsnp_GRCh38.ok'),
        'R/local/SNPlocs.Hsapiens.dbSNP{}.GRCh38/R/SNPlocs.Hsapiens.dbSNP{}.GRCh38.rdx'.format(dbsnp_v_GRCh38, dbsnp_v_GRCh38)
    log:
        'logs/install_bioconductor_dbsnp_GRCh38.log'
    shell:
        "("
        "mkdir -p R/local && "
        "R -e 'if (!requireNamespace(\"BiocManager\", quietly = TRUE)){{install.packages(\"BiocManager\", repos=\"https://cloud.r-project.org/\")}}; BiocManager::install(\"SNPlocs.Hsapiens.dbSNP{params[dbsnp_version]}.GRCh38\", lib=\"R/local/\");'"
        ") &> {log}"


# this parameter sets the version of the dbSNP R-package used by the rules below.
# as of 4.11.2021, this was the latest release available for GRCh37

dbsnp_v_GRCh37 = '144'

rule install_bioconductor_dbsnp_GRCh37:
    input:
        rules.install_bioconductor_dbsnp_GRCh38.output
    params:
        # version of dbSNP package to use, see https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP151.GRCh38.html
        dbsnp_version=dbsnp_v_GRCh37
    output:
        touch('R/install_dbsnp_GRCh37.ok'),
        'R/local/SNPlocs.Hsapiens.dbSNP{}.GRCh37/R/SNPlocs.Hsapiens.dbSNP{}.GRCh37.rdx'.format(dbsnp_v_GRCh37, dbsnp_v_GRCh37)
    log:
        'logs/install_bioconductor_dbsnp_GRCh37.log'
    shell:
        "("
        "mkdir -p R/local && "
        "R -e 'if (!requireNamespace(\"BiocManager\", quietly = TRUE)){{install.packages(\"BiocManager\", repos=\"https://cloud.r-project.org/\")}}; BiocManager::install(\"SNPlocs.Hsapiens.dbSNP{params[dbsnp_version]}.GRCh37\", lib=\"R/local/\");'"
        ") &> {log}"
        
rule intersect_1kg_hm3:
    input:
        bim_hm3=rules.make_hapmap3_bim.output,
        bim_1kg=rules.download_1kg.output['bim']
    output:
        mapping=temp("resources/1kg/1KGPhase_hm3_hg19_mapping_chr{chr}.tsv.gz")
    params:
        dbsnp_version=dbsnp_v_GRCh37,
        bim_hm3 = lambda wc, input: ','.join(input['bim_hm3'])
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    log:
        "logs/intersect_1kg_hm3/{chr}.log"
    shell:
        "("
        "{config[Rscript]} workflow/scripts/R/setup/intersect_1kg_hm3.R "
        "--bim {input[bim_1kg]} "
        "--hm3_bim {params[bim_hm3]} "
        "--chr {wildcards[chr]} "
        "--dbsnp_version {params[dbsnp_version]} "
        "--out {output[mapping]} "
        "--plink2 {config[plink2]} ;"
        ") &> {log} "
        
        
rule liftover_1kg_i_hm3_to_hg38:
    input:
        mapping=rules.intersect_1kg_hm3.output['mapping'],
        liftover=rules.download_liftover.output['liftover'],
        chain=rules.download_liftover.output['hg19_to_hg38_chain']
    output:
        bed_hg19=temp("resources/1kg/1KGPhase3_hm3_hg19_chr{chr}.BED"),
        bed_hg38=temp("resources/1kg/1KGPhase3_hm3_hg38_chr{chr}.BED"),
        unlifted="resources/1kg/1KGPhase3_hm3_hg38_chr{chr}.unlifted",
        liftover_hg38=temp("resources/1kg/1KGPhase3_hm3_hg38_liftover_chr{chr}.txt")
    log:
        "logs/liftover_1kg_i_hm3_to_hg38/{chr}.log"
    shell:
        "("
        "zcat {input[mapping]} | awk 'NR > 1' | awk -F '\t' 'BEGIN {{ chrom = 1; pos = 2; rsid_1kg = 3; rsid_dbsnp = 4; a1 = 5; a2 = 6;}} "
        "{{ print \"chr\"$chrom, $pos-1, $pos, $rsid_1kg, $rsid_dbsnp, $a1, $a2  }}' > {output.bed_hg19} && "
        "{input.liftover} {output.bed_hg19} {input.chain} {output.bed_hg38} {output.unlifted} && "
        "awk -F '\t' 'BEGIN {{OFS=\"\\t\"; chrom = 1; pos0 = 2; pos1 = 3; rsid_1kg = 4; rsid_dbsnp = 5; a1 = 6; a2 = 7; print \"chr\", \"pos_hg38\", \"rsid_1kg\", \"rsid_dbsnp\", \"a1\", \"a2\"}} "
        "{{ gsub(\"chr\", \"\", $chrom); print $chrom, $pos1, $rsid_1kg, $rsid_dbsnp, $a1, $a2 }}' {output.bed_hg38} > {output.liftover_hg38} "
        ") &> {log}"
        
        
rule generate_1kg_hm3_hg19_hg38_mapping:
    input:
        mapping=rules.intersect_1kg_hm3.output['mapping'],
        lifted=rules.liftover_1kg_i_hm3_to_hg38.output['liftover_hg38']
    output:
        mapping=temp('resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping_chr{chr}.tsv.gz')
    params:
        v_dbsnp_hg38=dbsnp_v_GRCh38,
        v_dbsnp_hg37=dbsnp_v_GRCh37
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    log:
        "logs/generate_1kg_hm3_hg19_hg38_mapping/{chr}.log"
    shell:
        "("
        "{config[Rscript]} workflow/scripts/R/setup/merge_1kg_hm3_hg19_hg38.R "
        "--mapping {input.mapping} "
        "--lifted {input.lifted} "
        "--chr {wildcards[chr]} "
        "--dbsnp_version_hg38 {params[v_dbsnp_hg38]} "
        "--dbsnp_version_hg37 {params[v_dbsnp_hg37]} "
        "--out {output[mapping]} "
        ")&>{log}"
        
        
rule all_generate_1kg_hm3_hg19_hg38_mapping:
    input:
        expand(rules.generate_1kg_hm3_hg19_hg38_mapping.output, chr=range(1,23))
        
        
##################################################
# filtering the 1000 Genomes genotypes - skipped #
##################################################

# Remo: these rules are here for documentation purposes only
#       these rules will be skipped, because we distibute their output with the pipeline.

rule filter_1kg_firstpass:
    # temporary first pass of filtering...
    # skipped.
    input:
        mapping=rules.generate_1kg_hm3_hg19_hg38_mapping.output['mapping'],
        bed=rules.download_1kg.output['bed'],
        fam=rules.download_1kg.output['fam'],
        bim=rules.download_1kg.output['bim']
    output:
        bed=temp('resources/1kg/tmp/1KGPhase3.chr{chr}.tmp.bed'),
        bim=temp('resources/1kg/tmp/1KGPhase3.chr{chr}.tmp.bim'),
        fam=temp('resources/1kg/tmp/1KGPhase3.chr{chr}.tmp.fam')
    params:
        out_prefix=lambda wc, output: output['bed'][:-4]
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    log:
        'logs/filter_1kg_firstpass/{chr}.log'
    shell:
        "{config[Rscript]} workflow/scripts/R/setup/position_based_hm3_harmonizer.R "
        "--bim {input[bim]} "
        "--out_prefix {params[out_prefix]} " 
        "--mapping {input[mapping]} "
        "--plink2 {config[plink2]} "
        "--rsid_col 'rsid_mrg' &> {log}"
        
        
rule allele_freq_1kg_pop_tmp:
    # 2.5 1000 Genomes allele frequency files
    # skipped.
    input:
        rules.create_ancestry.output,
        bed = expand(rules.filter_1kg_firstpass.output['bed'], chr=range(1,23)),
        bim = expand(rules.filter_1kg_firstpass.output['bim'], chr=range(1,23)),
        fam = expand(rules.filter_1kg_firstpass.output['fam'], chr=range(1,23))
    output:
        expand('resources/1kg/tmp/{{popul}}/chr{chr}.frq', chr=range(1,23), allow_missing=True)
    params:
        keep_file=lambda wc: "resources/1kg/keep_files/{}_samples.keep".format(wc['popul']),
        in_prefix=lambda wc, input: input['bed'][0].replace('1.tmp.bed',''),
        out_prefix=lambda wc, output: output[0][:-5]
    wildcard_constraints:
        pop='[A-Z]+'
    log:
        "logs/allele_freq_1kg_pop/{popul}.log"
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "for chrom in $(seq 1 22); do "
        "{config[plink1_9]} --bfile {params[in_prefix]}${{chrom}}.tmp "
        "--keep-allele-order "
        "--keep {params[keep_file]} "
        "--freq "
        "--out {params[out_prefix]}${{chrom}}; "
        "done "
        ") &> {log} "
        

rule run_allele_freq_1kg_superpop:
    # runs 2.5 for super populations
    # skipped.
    input:
        rules.create_ancestry.output,
        expand(rules.allele_freq_1kg_pop_tmp.output, popul=config['1kg_superpop'])


rule run_allele_freq_1kg_allancestry:
    # runs 2.5 for all ancestries combined
    # skipped.
    input:
        bed = expand(rules.filter_1kg_firstpass.output['bed'], chr=range(1,23)),
        bim = expand(rules.filter_1kg_firstpass.output['bim'], chr=range(1,23)),
        fam = expand(rules.filter_1kg_firstpass.output['fam'], chr=range(1,23))
    output:
        expand('resources/1kg/tmp/AllAncestry/chr{chr}.frq', chr=range(1,23))
    params:
        in_prefix=lambda wc, input: input['bed'][0].replace('1.tmp.bed',''),
        out_prefix=lambda wc, output: output[0][:-5]
    wildcard_constraints:
        pop='[A-Z]+'
    log:
        "logs/run_allele_freq_1kg_allancestry.log"
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "for chrom in $(seq 1 22); do "
        "{config[plink1_9]} --bfile {params[in_prefix]}${{chrom}}.tmp "
        "--keep-allele-order "
        "--freq "
        "--out {params[out_prefix]}${{chrom}}; "
        "done "
        ") &> {log} "
        
       
ruleorder: merge_1kg_hm3_mapping_with_maf_cached > merge_1kg_hm3_mapping_with_maf
 
rule merge_1kg_hm3_mapping_with_maf:
    # the maf threshold is applied here!
    # skipped.
    # the output of this rule is provided with the pipeline. This will skip the rules above marked with 'skipped'
    input:
        rules.run_allele_freq_1kg_allancestry.output,
        rules.run_allele_freq_1kg_superpop.input,
        expand(rules.generate_1kg_hm3_hg19_hg38_mapping.output, chr=range(1,23))
    output:
        'resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping.tsv.gz'
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    log:
        'logs/merge_1kg_hm3_mapping_with_maf.log'
    shell:
        "("
        "{config[Rscript]} workflow/scripts/R/setup/merge_1kg_hm3_mapping_with_maf.R "
        "--frq_path_prefix resources/1kg/tmp "
        "--mapping_prefix resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping_chr "
        "--out_prefix resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping "
        "--min_maf 0.01 "
        ") &> {log}"
        
rule merge_1kg_hm3_mapping_with_maf_cached:
    # skip the rules above and instead use the "cached" output
    input:
        ancient('resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping_cached.tsv.gz') 
    output:
        'resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping.tsv.gz'
    shell:
        'ln -s -r {input} {output}'

localrules: merge_1kg_hm3_mapping_with_maf_cached
        
#########################################
# filtering the 1000 Genomes genotypes  #
#########################################

# These are the rules that define the reference used in subsequent analyses

ruleorder: skip_1kg_hm3_processing > extract_hm3

rule skip_1kg_hm3_processing:
    output:
        expand("resources/1kg/1KGPhase3.w_hm3.chr{chrom}.{ext}", chrom=range(1,23), ext=['bed','bim','fam'])
    log:
        "logs/skip_1kg_hm3_processing.log"
    shell:
        "("
        "wget 'https://figshare.com/ndownloader/files/37059184?private_link=cea777bb772ed1dc4ca8' -O 1KGPhase3.w_hm3.chr.tar.gz && "
        "tar -xzvf 1KGPhase3.w_hm3.chr.tar.gz && rm 1KGPhase3.w_hm3.chr.tar.gz && "
        "find resources/1kg/ -name '1KGPhase3.w_hm3.chr*' -type f -exec touch {{}} + "
        ") &> {log} "

rule extract_hm3:
    # Outputs the 1000 Genomes filtered to HapMap3 variants
    # Generates per-chromosome files.
    input:
        mapping='resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping.tsv.gz',
        bed=rules.download_1kg.output['bed'],
        fam=rules.download_1kg.output['fam'],
        bim=rules.download_1kg.output['bim']
    output:
        bed='resources/1kg/1KGPhase3.w_hm3.chr{chr}.bed',
        bim='resources/1kg/1KGPhase3.w_hm3.chr{chr}.bim',
        fam='resources/1kg/1KGPhase3.w_hm3.chr{chr}.fam'
    params:
        out_prefix=lambda wc, output: output['bed'][:-4]
    resources:
        time="03:00:00",
        mem_mb=8000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    log:
        'logs/extract_hm3/{chr}.log'
    shell:
        "{config[Rscript]} workflow/scripts/R/setup/position_based_hm3_harmonizer.R "
        "--bim {input[bim]} "
        "--out_prefix {params[out_prefix]} " 
        "--mapping {input[mapping]} "
        "--plink2 {config[plink2]} "
        "--tmpdir ./temp "
        "--rsid_col 'rsid' &> {log}"
    
    
rule extract_hm3_gw:
    # Outputs the 1000 Genomes filtered to HapMap3 variants
    # Generates a single file with all chromosomes
    input:
        expand(rules.extract_hm3.output, chr=range(1,23))
    output:
        bed='resources/1kg/1KGPhase3.w_hm3.GW.bed',
        bim='resources/1kg/1KGPhase3.w_hm3.GW.bim',
        fam='resources/1kg/1KGPhase3.w_hm3.GW.fam'
    log:
        "logs/extract_hm3_gw.log"
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    shell:
        "("
        "ls resources/1kg/1KGPhase3.w_hm3.chr*.bed | sed -e 's/\.bed//g' > resources/1kg/merge_list.txt ;"
        "{config[plink1_9]} --merge-list resources/1kg/merge_list.txt "
        "--make-bed "
        "--out resources/1kg/1KGPhase3.w_hm3.GW "
        ") &> {log} "
        

############################################
# Allele frequency calculation - plink 1.9 #
############################################

rule allele_freq_pop:
    # 2.5 1000 Genomes allele frequency files
    # Here we create files containing ancestry specific minor allele frequency estimates for the HapMap3 SNPs based on the 1000 Genomes Phase 3 data. This information is mainly used for mean-imputation of missing SNPs during genotype-based scoring. This avoids target sample specific minor allele frequencies being used for mean imputation which may not be available, and will vary between target samples.
    input:
        rules.create_ancestry.output,
        expand(rules.extract_hm3.output, chr=range(1,23), allow_missing=True)
    output:
        pop_list=expand("resources/1kg/freq_files/{{popul}}/1KGPhase3.w_hm3.{{popul}}.chr{chr}.frq", chr=range(1,23))
    wildcard_constraints:
        pop='[A-Z]+'
    log:
        "logs/allele_freq_pop/{popul}.log"
    singularity:
        config['singularity']['all']
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    shell:
        "("
        "for chrom in $(seq 1 22); do "
        "{config[plink1_9]} --bfile resources/1kg/1KGPhase3.w_hm3.chr${{chrom}} "
        "--keep resources/1kg/keep_files/{wildcards[popul]}_samples.keep "
        "--freq "
        "--out resources/1kg/freq_files/{wildcards[popul]}/1KGPhase3.w_hm3.{wildcards[popul]}.chr${{chrom}}; "
        "done "
        ") &> {log} "


rule run_allele_freq_pop:
    # runs 2.5 for populations
    input:
        rules.create_ancestry.output,
        expand(rules.allele_freq_pop.output, popul=config['1kg_pop'])
    output:
        touch("resources/1kg/freq_files/all_pop.ok")
        
        
rule run_allele_freq_superpop:
    input:
        rules.create_ancestry.output,
        expand(rules.allele_freq_pop.output, popul=config['1kg_superpop'])
    output:
        touch("resources/1kg/freq_files/all_superpop.ok")
        

rule run_allele_freq_allancestry:
    # runs 2.5 for all ancestries combined
    input:
        expand(rules.extract_hm3.output, chr=range(1,23), allow_missing=True)
    output:
        expand("resources/1kg/freq_files/AllAncestry/1KGPhase3.w_hm3.AllAncestry.chr{chr}.frq", chr=range(1, 23))
    log:
        "logs/run_allele_freq_allancestry.log"
    singularity:
        config['singularity']['all']
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    shell:
        "("
        "for chr in $(seq 1 22); do "
        "{config[plink1_9]} --bfile resources/1kg/1KGPhase3.w_hm3.chr${{chr}} "
        "--freq "
        "--out resources/1kg/freq_files/AllAncestry/1KGPhase3.w_hm3.AllAncestry.chr${{chr}}; "
        "done "
        ") &> {log} "


rule run_allele_freq_all:
    # this rule executes all the rules above
    input:
        rules.run_allele_freq_superpop.output,
        rules.run_allele_freq_allancestry.output
    output:
        touch("resources/1kg/freq_files/all_plink1.ok")


##########################################
# Allele frequency calculation - plink 2 #
##########################################

# Plink 1.x allele frequency files lead to segmentation faults when combined with plink2 scoring (?)
# For this reason we create plink2 .afreq files to use for polygenic scoring.

rule allele_freq_pop_plink2:
    input:
        rules.create_ancestry.output,
        expand(rules.extract_hm3.output, chr=range(1,23), allow_missing=True)
    output:
        pop_list=expand("resources/1kg/freq_files/{{popul}}/1KGPhase3.w_hm3.{{popul}}.chr{chr}.afreq", chr=range(1,23))
    wildcard_constraints:
        pop='[A-Z]+'
    log:
        "logs/allele_freq_pop_plink2/{popul}.log"
    resources:
        time="01:00:00",
        mem_mb=4000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "for chrom in {{1..22}}; do "
        "{config[plink2]} --bfile resources/1kg/1KGPhase3.w_hm3.chr${{chrom}} "
        "--keep resources/1kg/keep_files/{wildcards[popul]}_samples.keep "
        "--freq "
        "--out resources/1kg/freq_files/{wildcards[popul]}/1KGPhase3.w_hm3.{wildcards[popul]}.chr${{chrom}}; "
        "done "
        ") &> {log} "


rule run_allele_freq_pop_plink2:
    input:
        rules.create_ancestry.output,
        expand(rules.allele_freq_pop_plink2.output, popul=config['1kg_pop'])
    output:
        touch("resources/1kg/freq_files/all_pop_plink2.ok")

rule run_allele_freq_superpop_plink2:
    input:
        rules.create_ancestry.output,
        expand(rules.allele_freq_pop_plink2.output, popul=config['1kg_superpop'])
    output:
        touch("resources/1kg/freq_files/all_superpop_plink2.ok")


rule run_allele_freq_allancestry_plink2:
    input:
        expand(rules.extract_hm3.output, chr=range(1,23), allow_missing=True)
    output:
        expand("resources/1kg/freq_files/AllAncestry/1KGPhase3.w_hm3.AllAncestry.chr{chr}.afreq", chr=range(1, 23))
    log:
        "logs/run_allele_freq_allancestry.log"
    resources:
        time="01:00:00",
        mem_mb=4000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "for chr in {{1..22}}; do "
        "{config[plink2]} --bfile resources/1kg/1KGPhase3.w_hm3.chr${{chr}} "
        "--freq "
        "--out resources/1kg/freq_files/AllAncestry/1KGPhase3.w_hm3.AllAncestry.chr${{chr}}; "
        "done "
        ") &> {log} "


rule run_allele_freq_all_plink2:
    # this rule executes all the rules above
    input:
        rules.run_allele_freq_superpop_plink2.output,
        rules.run_allele_freq_allancestry_plink2.output
    output:
        touch("resources/1kg/freq_files/all.ok")


##########################
# Ancestry scoring - PCA #
##########################

# create files needed for ancestry scoring

rule ancestry_scoring:
    # 3 Ancestry scoring
    # https://opain.github.io/GenoPred/Pipeline_prep.html
    # In this section we perform principal components analysis of genotypic data in the 1000 Genomes reference to identify the main axes of population structure, which typically correspond to ancestral differences. We calculate these principal components (PCs) of ancestry across the full reference, and within super populations to detect broad and ancestry-specific axes of variance. After PCA, we idenitfy which variants are associated with the PCs (SNP-weights) to calculate ancestry scores on the same axes of variance in future target samples. This can be used to infer the ancestry of an individual which is an important factor to consider when performing genotype-baed prediction.
    # The script will always create .eigenvec and .eigenvec.var files. These contain the PCs score for individuals in the reference dataset, and the SNP-weights for the PCs respectively.
    input:
        rules.create_ancestry.output,
        expand(rules.extract_hm3.output, chr=range(1,23), allow_missing=True)
    output:
        eigenvec="resources/1kg/Score_files_for_ancestry/{popul}/1KGPhase3.w_hm3.{popul}.eigenvec",
        pc_scale="resources/1kg/Score_files_for_ancestry/{popul}/1KGPhase3.w_hm3.{popul}.scale",
        eigenvec_var="resources/1kg/Score_files_for_ancestry/{popul}/1KGPhase3.w_hm3.{popul}.eigenvec.var"
    log:
        "logs/ancestry_scoring/{popul}.log"
    wildcard_constraints:
        popul='[A-Z]+'
    resources:
        time="08:00:00",
        mem_mb=16000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "{config[Rscript]} {config[GenoPred_dir]}/Scripts/ancestry_score_file_creator/ancestry_score_file_creator.R "
        "--ref_plink_chr resources/1kg/1KGPhase3.w_hm3.chr "
        "--ref_keep resources/1kg/keep_files/{wildcards[popul]}_samples.keep "
        "--plink {config[plink1_9]} "
        "--plink2 {config[plink2]} "
        "--n_pcs 100 "
        "--output resources/1kg/Score_files_for_ancestry/{wildcards[popul]}/1KGPhase3.w_hm3.{wildcards[popul]}"
        ") &> {log}"


rule ancestry_scoring_allancestry:
    # 3 Ancestry scoring
    # https://opain.github.io/GenoPred/Pipeline_prep.html
    # If --ref_pop_scale is specified, the script will also creates files stating the mean and standard deviation of the PCs for each group. Furthermore, it will derive an elastic net model predicting each group, and report the accuracy of the derived models.
    input:
        rules.create_ancestry.output,
        expand(rules.extract_hm3.output, chr=range(1,23), allow_missing=True)
    output:
        pop_enet_model="resources/1kg/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.pop_enet_model.rds",
        eigenvec="resources/1kg/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.eigenvec",
        pc_scale="resources/1kg/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.scale",
        eigenvec_var="resources/1kg/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.eigenvec.var"
    log:
        'logs/ancestry_scoring_allancestry.log'
    resources:
        time="08:00:00",
        mem_mb=32000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "{config[Rscript]} {config[GenoPred_dir]}/Scripts/ancestry_score_file_creator/ancestry_score_file_creator.R "
        "--ref_plink_chr resources/1kg/1KGPhase3.w_hm3.chr "
        "--plink {config[plink1_9]} "
        "--plink2 {config[plink2]} "
        "--n_pcs 100 "
        "--ref_pop_scale resources/1kg/super_pop_keep.list "
        "--output resources/1kg/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry"
        ") &> {log}"

