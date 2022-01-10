# Remo:
# rules related to genotype harmonization including pre-processing of HapMap3 variants and intersection with 1000 Genomes.
# this module should eventually replace the standard pre-processing based on the the ldsc-hapmap3 list (which uses rsIDs only...)

######################################
# Download the 1000Genomes reference #
######################################

rule download_1kg:
    # generate 1000 Genomes PLINK files
    # Download vcf file, convert to plink
    output:
        bim="resources/1kg/1KGPhase3.chr{chr}.bim",
        bed="resources/1kg/1KGPhase3.chr{chr}.bed",
        fam="resources/1kg/1KGPhase3.chr{chr}.fam"
    params:
        # v5b - no rsIDs!: ftp_path=lambda wc: "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz".format(wc['chr']),
        # v5 - has rsIDs:
        ftp_path=lambda wc: "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140708_previous_phase3/v5_vcfs/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz".format(wc['chr']),
        out = lambda wc, output: output['bed'][:-4]
    log:
        "logs/download_1kg/{chr}.log"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "wget {params[ftp_path]} -nv -O resources/1kg/chr{wildcards[chr]}.vcf.gz && "
        "gunzip resources/1kg/chr{wildcards[chr]}.vcf.gz && "
        "{config[plink1_9]} --vcf resources/1kg/chr{wildcards[chr]}.vcf "
        "--make-bed "
        "--out {params[out]} && "
        "rm resources/1kg/chr{wildcards[chr]}.vcf && "
        "rm {params[out]}.log "
        ") &> {log} "


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
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    log:
        "logs/intersect_1kg_hm3/{chr}.log"
    shell:
        "("
        "Rscript workflow/scripts/R/setup/intersect_1kg_hm3.R "
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
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    log:
        "logs/generate_1kg_hm3_hg19_hg38_mapping/{chr}.log"
    shell:
        "("
        "Rscript workflow/scripts/R/setup/merge_1kg_hm3_hg19_hg38.R "
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
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    log:
        'logs/filter_1kg_firstpass/{chr}.log'
    shell:
        "Rscript workflow/scripts/R/setup/position_based_hm3_harmonizer.R "
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
        keep_file=lambda wc: "{}/keep_files/{}_samples.keep".format(config['Geno_1KG_dir'], wc['popul']),
        in_prefix=lambda wc, input: input['bed'][0].replace('1.tmp.bed',''),
        out_prefix=lambda wc, output: output[0][:-5]
    wildcard_constraints:
        pop='[A-Z]+'
    log:
        "logs/allele_freq_1kg_pop/{popul}.log"
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
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    log:
        'logs/merge_1kg_hm3_mapping_with_maf.log'
    shell:
        "("
        "Rscript workflow/scripts/R/setup/merge_1kg_hm3_mapping_with_maf.R "
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

        
#########################################
# filtering the 1000 Genomes genotypes  #
#########################################
        
rule extract_hm3:
    # the second and final filtering pass
    # replaces the output of the original "extract_hm3"
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
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    log:
        'logs/filter_1kg_secondpass/{chr}.log'
    shell:
        "Rscript workflow/scripts/R/setup/position_based_hm3_harmonizer.R "
        "--bim {input[bim]} "
        "--out_prefix {params[out_prefix]} " 
        "--mapping {input[mapping]} "
        "--plink2 {config[plink2]} "
        "--rsid_col 'rsid' &> {log}"
    
    
rule extract_hm3_gw:
    input:
        expand(rules.extract_hm3.output, chr=range(1,23))
    output:
        bed='resources/1kg/1KGPhase3.w_hm3.GW.bed',
        bim='resources/1kg/1KGPhase3.w_hm3.GW.bim',
        fam='resources/1kg/1KGPhase3.w_hm3.GW.fam'
    log:
        "logs/generate_gw_1kg_hm3_ref.log"
    shell:
        "("
        "ls resources/1kg/1KGPhase3.w_hm3.chr*.bed | sed -e 's/\.bed//g' > resources/1kg/merge_list.txt ;"
        "{config[plink1_9]} --merge-list resources/1kg/merge_list.txt "
        "--make-bed "
        "--out resources/1kg/1KGPhase3.w_hm3.GW "
        ") &> {log} "
        
        
rule allele_freq_pop:
    # 2.5 1000 Genomes allele frequency files
    # Here we create files containing ancestry specific minor allele frequency estimates for the HapMap3 SNPs based on the 1000 Genomes Phase 3 data. This information is mainly used for mean-imputation of missing SNPs during genotype-based scoring. This avoids target sample specific minor allele frequencies being used for mean imputation which may not be available, and will vary between target samples.
    input:
        rules.create_ancestry.output,
        expand(rules.extract_hm3.output, chr=range(1,23), allow_missing=True)
    params:
        keep_file=lambda wc: "{}/keep_files/{}_samples.keep".format(config['Geno_1KG_dir'], wc['popul'])
    output:
        pop_list=expand("{}/freq_files/{{popul}}/1KGPhase3.w_hm3.{{popul}}.chr{{chr}}.frq".format(config['Geno_1KG_dir']), chr=range(1,23), allow_missing=True)
    wildcard_constraints:
        pop='[A-Z]+'
    log:
        "logs/allele_freq_pop/{popul}.log"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "for chrom in $(seq 1 22); do "
        "{config[plink1_9]} --bfile {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.chr${{chrom}} "
        "--keep {config[Geno_1KG_dir]}/keep_files/{wildcards[popul]}_samples.keep "
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
        expand(rules.extract_hm3.output, chr=range(1,23), allow_missing=True)
    output:
        expand("{}/freq_files/AllAncestry/1KGPhase3.w_hm3.AllAncestry.chr{{chr}}.frq".format(config['Geno_1KG_dir']), chr=range(1, 23))
    log:
        "logs/run_allele_freq_allancestry/run_allele_freq_allancestry.log"
    singularity:
        config['singularity']['all']
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
    output:
        touch("{}/freq_files/all.ok".format(config['Geno_1KG_dir']))


rule ancestry_scoring:
    # 3 Ancestry scoring
    # https://opain.github.io/GenoPred/Pipeline_prep.html
    # In this section we perform principal components analysis of genotypic data in the 1000 Genomes reference to identify the main axes of population structure, which typically correspond to ancestral differences. We calculate these principal components (PCs) of ancestry across the full reference, and within super populations to detect broad and ancestry-specific axes of variance. After PCA, we idenitfy which variants are associated with the PCs (SNP-weights) to calculate ancestry scores on the same axes of variance in future target samples. This can be used to infer the ancestry of an individual which is an important factor to consider when performing genotype-baed prediction.
    # The script will always create .eigenvec and .eigenvec.var files. These contain the PCs score for individuals in the reference dataset, and the SNP-weights for the PCs respectively.
    input:
        rules.create_ancestry.output,
        expand(rules.extract_hm3.output, chr=range(1,23), allow_missing=True)
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
    resources:
        time="08:00:00",
        mem_mb=16000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    singularity:
        config['singularity']['all']
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
        expand(rules.extract_hm3.output, chr=range(1,23), allow_missing=True)
    output:
        pop_enet_model="{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.pop_enet_model.rds".format(config['Geno_1KG_dir']),
        eigenvec="{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.eigenvec".format(config['Geno_1KG_dir']),
        pc_scale="{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.scale".format(config['Geno_1KG_dir']),
        eigenvec_var="{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.eigenvec.var".format(config['Geno_1KG_dir'])
    log:
        # the script is configured to write log files here
        "{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.log".format(config['Geno_1KG_dir'])
    resources:
        time="08:00:00",
        mem_mb=32000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    singularity:
        config['singularity']['all']
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
    # the Polygenic scoring steps (>=4) that follow are implemented in different modules
    input:
        '{}/Score_files_for_ancestry/EUR/1KGPhase3.w_hm3.EUR.eigenvec.var'.format(config['Geno_1KG_dir']),
        '{}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.eigenvec.var'.format(config['Geno_1KG_dir']),
        rules.run_allele_freq_all.output,
        rules.download_liftover.output
        
        
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
        geno = lambda wc: get_input_file_from_type(target_list.loc[wc['bbid'], 'path'],wc['chr'],target_list.loc[wc['bbid'], 'type']),
        liftover=rules.download_liftover.output['liftover'],
        chain=rules.download_liftover.output['hg19_to_hg38_chain'],
        plink2=config['plink2'],
        plink19=config['plink1_9']
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
        mem_mb=8000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    log:
        'logs/harmonize_target_genotypes/{bbid}/{chr}.log'
    shell:
        '('
        'Rscript workflow/scripts/GenoPred/Scripts/geno_to_plink/geno_to_plink.R '
        '--target {params[in_prefix_target]} '
        '--ref {params[in_prefix_ref]} '
        '--format {params[arg][0]} '
        '--plink {config[plink1_9]} '
        '--plink2 {config[plink2]} '
        '--qctool2 {config[qctool2]} '
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
