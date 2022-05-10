
rule install_basics:
    # download plink 1.9, plink 2, and clone GenoPred and ldsc
    output:
        "bin/plink",
        "bin/plink2",
        directory('workflow/scripts/GenoPred'),
        directory('workflow/scripts/ldsc'),
        "bin/qctool"
    shell:
        "bash install_basics.sh"
        

rule download_hapmap3_snplist:
    input:
        'resources/1kg/1KGPhase3_hm3_hg19_hg38_mapping.tsv.gz'
    output:
        "resources/HapMap3_snplist/w_hm3.snplist"
    log:
        "logs/download_hapmap3_snplist.log"
    shell:
        "("
        "zcat {input} | cut -f2-4 | "
        ' awk \'BEGIN{{OFS="\t"; print "SNP", "A1", "A2"}}{{if(NR > 1){{print $0}}}}\' > {output} '
        ") &> {log}"
        
        # old version: use the LDSC list
        # "cd resources/HapMap3_snplist && "
        # "wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 && "
        # "bunzip2 w_hm3.snplist.bz2"


rule initialize_synthetic:
    # unzip the synthetic data
    # rename the columns in the sumstats file
    output:
        bed=expand('resources/synthetic_data/synthetic250_chr{chr}.bed', chr=range(1,23)),
        bim=expand('resources/synthetic_data/synthetic250_chr{chr}.bim', chr=range(1,23)),
        fam=expand('resources/synthetic_data/synthetic250_chr{chr}.fam', chr=range(1,23)),
        ss='resources/synthetic_data/Herr0.3_pPoly0.005.PHENO1.glm.linear.renamed.gz',
        pheno='custom_input/synth/phenotypes/synthetic01.tsv.gz'
    shell:
        "(cd resources/synthetic_data/ && unzip synthetic_genotypes_250.zip && "
        "zcat Herr0.3_pPoly0.005.PHENO1.glm.linear.gz | awk 'BEGIN{{OFS=\"\t\"; print \"CHR\", \"POS\", \"SNP\", \"A2\", \"A1\", \"N\", \"BETA\", \"SE\", \"P\"}}{{if(NR>1){{print $1, $2, $3, $4, $5, $8, $9, $10, $14}}}}' | gzip > Herr0.3_pPoly0.005.PHENO1.glm.linear.renamed.gz ); "
        "ln -s -r \"$(readlink -f resources/synthetic_data/pheno250.tsv.gz)\" {output[pheno]}"


rule download_test_data:
    # rule to download some small test data from figshare (simple PRS for 3 methods and 10 phenotypes)
    # touches the output files to fool snakemake's timestamp checks
    # note: this is incompatible with the rules below ("download_prs_for_methods_coparison_may2022", "unpack_prs_for_methods_coparison_may2022"). (outputs will be overwritten)
    # meant for use with config/studies_for_methods_comparison.tsv
    shell:
        "( mkdir -p resources/test_prs/ && "
        "wget -O resources/test_prs/PRS.tar.gz https://figshare.com/ndownloader/files/33905531?private_link=8ac0f09450555d6ba6dd && "
        "cd resources/test_prs/ && "
        "tar -xvf PRS.tar.gz && "
        "find . -type f -exec touch {{}} + ); " # touching the extracted files so the snakemake timestamp checks work
        "mkdir -p prs && cd prs; "
        "for DIR in ../resources/test_prs/*/; do "
        "METHOD=${{DIR##../resources/test_prs/}}; "
        "mkdir -p $METHOD; cd $METHOD; "
        "for SUBDIR in ../${{DIR}}/*; "
        "do ln -s -r $SUBDIR; done;" # linking the directories to the standard directory for polygenic scores
        "cd ../; done "


rule download_prs_for_methods_coparison_may2022:
    # rule to download pre-computed PRS scoring files
    # note: this is incompatible with the rule above. (outputs will be overwritten)
    output:
        'prs.tar.gz'
    shell:
        "wget -O prs.tar.gz https://figshare.com/ndownloader/files/35016199?private_link=7f738b939e1cba580708"



with open('config/studies_for_methods_comparison.tsv', 'r') as infile:
    may2022_study_ids = [x.split('\t')[0] for x in infile][1:]


rule unpack_prs_for_methods_coparison_may2022:
    # rule to unpack pre-computed PRS scoring files
    # note: this is incompatible with the "download_test_data" rule above. (outputs will be overwritten)
    # touches the output files to fool snakemake's timestamp checks
    # meant for use with config/studies_for_methods_comparison.tsv
    input:
        'prs.tar.gz'
    params:
        study_ids = ' '.join(may2022_study_ids),
        prs_methods = ' '.join(config['prs_methods'])
    shell:
        'tar -xvf prs.tar.gz; '
        'for study in {params[study_ids]}; do '
        'for method in {params[prs_methods]}; do '
        'touch "prs/${{method}}/${{study}}/ok" && '
        'find "prs/${{method}}/${{study}}" -type f -exec touch {{}} +; '
        'done; '
        'done'


localrules:
    download_prs_for_methods_coparison_may2022,
    unpack_prs_for_methods_coparison_may2022


rule cleanup_after_setup:
    # removes unnecessary intermediate output files
    log:
        "logs/cleanup.log"
    shell:
        "("
        "rm -f resources/1kg/ALL.chr*.phase3_shapeit2_mvncall_integrated_v*.*.genotypes.vcf.gz ; "
        "rm -f resources/1kg/merge_list.txt ; "
        "rm -f resources/1kg/1KGPhase3.chr*.bed ; "
        "rm -f resources/1kg/1KGPhase3.chr*.bim ; "
        "rm -f resources/1kg/1KGPhase3.chr*.fam ; "
        "rm -f resources/1kg/1KGPhase3.chr*.nosex ; "
        "rm -f resources/1kg/1KGPhase3.chr*.log ; "
        "rm -f resources/1kg/1KGPhase3.chr*.extract ; "
        "rm -f resources/1kg/1KGPhase3.w_hm3.chr*.nosex ; "
        "rm -f resources/1kg/1KGPhase3.w_hm3.chr*.log ; "
        "rm -f resources/1kg/freq_files/*/*log ; "
        "rm -f resources/1kg/freq_files/*/*noseq ; "
        ") &> {log}"
