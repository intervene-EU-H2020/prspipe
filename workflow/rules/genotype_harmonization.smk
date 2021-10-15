
# Remo:
# rules related to genotype harmonization including pre-processing of HapMap3 variants and intersection with 1000 Genomes.
# this module should eventually replace the standard pre-processing based on the the ldsc-hapmap3 list (which uses rsIDs only...)

wildcard_constraints:
    chr="[0-9]+"

###########################################################
# This first part serves only documentation purposes.     #
# The outputs of this part already ship with the pipeline #   
###########################################################


rule download_liftover:
    # Note: won't be necessary if run with singularity
    output:
        liftover="liftover/liftOver",
        hg18_to_hg19_chain="liftover/hg18ToHg19.over.chain",
        hg18_to_hg38_chain="liftover/hg18ToHg38.over.chain"
    log:
        "logs/download_liftover.log"
    shell:
        "("
        "mkdir -p liftover && "
        "cd liftover && "
        "wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && "
        "wget https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz && "
        "gunzip hg18ToHg19.over.chain.gz && "
        "wget https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz && "
        "gunzip hg18ToHg38.over.chain.gz && "
        "chmod u+x liftOver "
        ") &> {log} "


rule download_hapmap3:
    output:
        ped="resources/hapmap3/hapmap3_r3_b36_fwd.consensus.qc.poly.ped",
        map="resources/hapmap3/hapmap3_r3_b36_fwd.consensus.qc.poly.map"
    log:
        "logs/download_hapmap3.log"
    shell:
        "("
        "mkdir -p resources/hapmap3 && "
        "cd resources/hapmap3 && "
        "wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.ped.gz && "
        "wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.map.gz && "
        "gunzip *.gz "
        ") &> {log}"
        
        
rule make_hapmap3_bim:
    input:
        rules.download_hapmap3.output
    output:
        bim="resources/hapmap3/hapmap3.bim"
    log:
        "logs/make_hapmap3_bim.log"
    shell:
        "("
        "{config[plink1_9]} --ped {input.ped} "
        "--map {input.map} "
        "--out resources/hapmap3/hapmap3 "
        "--make-just-bim "
        ") &> {log}"
        
        
rule liftover_hapmap3_to_hg19:
    input:
        bim=rules.make_hapmap3_bim.output['bim'],
        liftover=rules.download_liftover.output['liftover'],
        chain=rules.download_liftover.output['hg18_to_hg19_chain']
    output:
        bed_hg18=temp("resources/hapmap3/hapmap3.BED"),
        bed_hg19=temp("resources/hapmap3/hapmap3_hg19.BED"),
        unlifted="resources/hapmap3/hapmap3_hg19.unlifted",
        hapmap_hg19="resources/hapmap3/hapmap3_hg19_liftover.txt"
    log:
        "logs/liftover_hapmap3_to_hg19.log"
    shell:
        "("
        "awk -F '\t' 'BEGIN {{ chrom = 1; rsid = 2; dist = 3; pos = 4; a1 = 5; a2 = 6;}} "
        "{{ print \"chr\"$chrom, $pos-1, $pos, $rsid, $a1, $a2  }}' {input.bim} > {output.bed_hg18} && "
        "{input.liftover} {output.bed_hg18} {input.chain} {output.bed_hg19} {output.unlifted} && "
        "awk -F '\t' 'BEGIN {{ chrom = 1; pos0 = 2; pos1 = 3; rsid = 4; a1 = 5; a2 = 6 ; print \"chr\", \"pos_hg19\", \"rsid\", \"a1\", \"a2\"}} "
        "{{ gsub(\"chr\", \"\", $chrom); print $chrom, $pos1, $rsid, $a1, $a2 }}' {output.bed_hg19} > {output.hapmap_hg19} "
        ") &> {log}"
        
    
rule liftover_hapmap3_to_hg38:
    input:
        bim=rules.make_hapmap3_bim.output['bim'],
        liftover=rules.download_liftover.output['liftover'],
        chain=rules.download_liftover.output['hg18_to_hg38_chain']
    output:
        bed_hg18=temp("resources/hapmap3/hapmap3.BED"),
        bed_hg38=temp("resources/hapmap3/hapmap3_hg38.BED"),
        unlifted="resources/hapmap3/hapmap3_hg38.unlifted",
        hapmap_hg38="resources/hapmap3/hapmap3_hg38_liftover.txt"
    log:
        "logs/liftover_hapmap3_to_hg38.log"
    shell:
        "("
        "awk -F '\t' 'BEGIN {{ chrom = 1; rsid = 2; dist = 3; pos = 4; a1 = 5; a2 = 6;}} "
        "{{ print \"chr\"$chrom, $pos-1, $pos, $rsid, $a1, $a2  }}' {input.bim} > {output.bed_hg18} && "
        "{input.liftover} {output.bed_hg18} {input.chain} {output.bed_hg38} {output.unlifted} && "
        "awk -F '\t' 'BEGIN {{ chrom = 1; pos0 = 2; pos1 = 3; rsid = 4; a1 = 5; a2 = 6 ; print \"chr\", \"pos_hg38\", \"rsid\", \"a1\", \"a2\"}} "
        "{{ gsub(\"chr\", \"\", $chrom); print $chrom, $pos1, $rsid, $a1, $a2 }}' {output.bed_hg38} > {output.hapmap_hg38} "
        ") &> {log}"


rule merge_hapmap3_liftover_output:
    # requires R data.table
    # the output of this step is stored, so the previous steps can be skipped.
    input:
        hmhg19=rules.liftover_hapmap3_to_hg19.output['hapmap_hg19'],
        hmhg38=rules.liftover_hapmap3_to_hg38.output['hapmap_hg38'],
        hmhg18_bim=rules.make_hapmap3_bim.output
    output:
        "resources/hapmap3/hapmap3_mapping.tsv.gz"
    log:
        "logs/merge_hapmap3_liftover_output.log"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "Rscript workflow/scripts/R/setup/merge_hapmap3_positions.R "
        "{input.hmhg18_bim} "
        "{hmhg19} "
        "{hmhg38} "
        "{output} "
        ") &> {log}"



#######################################################
# The part below should be run in the single biobanks #
#######################################################


rule download_1kg:
    # generate 1000 Genomes PLINK files
    # Download vcf file, convert to plink
    output:
         bed="resources/1kg/1KGPhase3.chr{{chr}}.bed".format(config['Geno_1KG_dir']),
         fam="resources/1kg/1KGPhase3.chr{{chr}}.fam".format(config['Geno_1KG_dir']),
         bim="resources/1kg/1KGPhase3.chr{{chr}}.bim".format(config['Geno_1KG_dir'])
    params:
        ftp_path=lambda wc: "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz".format(wc['chr'])
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
        "--out resources/1kg/1KGPhase3.chr{wildcards[chr]} && "
        "rm resources/1kg/chr{wildcards[chr]}.vcf"
        ") &> {log} "
        
        
rule intersect_1kg_hm3:
    # requires R data.table and optparse + plink2
    input:
        bim=rules.download_1kg.output['bim'],
        mapping=rules.merge_hapmap3_liftover_output.output
    output:
        bed='resources/1kg/1KGPhase3.chr{chr}.w_hm3.bed',
        bim='resources/1kg/1KGPhase3.chr{chr}.w_hm3.bim',
        fam='resources/1kg/1KGPhase3.chr{chr}.w_hm3.fam'
    params:
        # output prefix
        out=lambda wc, output: output['bed'][:-4] 
    log:
        "logs/intersect_1kg_hm3/{chr}.log"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "Rscript workflow/scripts/R/setup/position_based_hm3_harmonizer.R "
        "--bim {input.bim} "
        "--genome_build hg19 "
        "--mapping {input.mapping} "
        "--plink2 {config[plink2]} "
        "--out_prefix {params[out]} "
        ") &> {log}"

rule all_intersect_1kg_hm3:
    input:
        expand(rules.intersect_1kg_hm3.output, chr=range(1,23))
        
        
rule merge_1kg_hm3_gw:
    input:
        expand(rules.intersect_1kg_hm3.output, chr=range(1,23))
    output:
        bed='resources/1kg/1KGPhase3.GW.w_hm3.bed',
        bim='resources/1kg/1KGPhase3.GW.w_hm3.bim',
        fam='resources/1kg/1KGPhase3.GW.w_hm3.fam'
    log:
        'logs/merge_1kg_hm3_gw.log'
    shell:
        "("
        "ls resources/1kg/1KGPhase3.w_hm3.chr*.bed | sed -e 's/\.bed//g' > resources/1kg/merge_list.txt ;"
        "{config[plink2]} --merge-list resources/1kg/merge_list.txt "
        "--make-bed "
        "--out {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        ") &> {log}"
    
    
