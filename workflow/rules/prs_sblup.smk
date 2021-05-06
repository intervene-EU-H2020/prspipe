# Rules for preparing score and scale files for polygenic scoring using sblup

# TODO add automatic download for ld_ref - currently done manually
# TODO currently uses eur_w_ld_chr as LD reference panel downloaded from https://alkesgroup.broadinstitute.org/LDSCORE/ - need to replace with our own? or support different ancestries at least?

rule sblup_prep:
    # Note that LDSC requires python 2 so Snakemake will setup this environment using the given .yml file
    input: 
        "{}/{{study}}.{{ancestry}}.cleaned.gz".format(config['Base_sumstats_dir'])
    output:
        "{}/sblup/{{study}}/1KGPhase3.w_hm3.{{ancestry}}.{{study}}".format(config['Base_sumstats_dir'])
    log:
        "logs/base_sumstats/prs_sblup_{study}.{ancestry}.log"
    conda:
        "../../{}/environment.yml".format(config['LDSC_dir'])
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/polygenic_score_file_creator_SBLUP/polygenic_score_file_creator_SBLUP.R "
        "--ref_plink {config[Geno_1KG_dir]}/1KGPhase3.w_hm3.GW "
        "--ref_keep {config[Geno_1KG_dir]}/keep_files/{wildcards[ancestry]}_samples.keep "
        "--sumstats {input} "
        "--plink {config[plink1_9]} "
        "--gcta {config[gcta]} "
        "--munge_sumstats {config[LDSC_dir]}/munge_sumstats.py "
        "--ldsc {config[LDSC_dir]}/ldsc.py "
        "--ldsc_ref {config[ldsc_ref]} "
        "--hm3_snplist {config[HapMap3_snplist_dir]}/w_hm3.snplist "
        "--memory 50000 "
        "--n_cores 6 "
        "--output {output} "
        "--ref_pop_scale {config[Geno_1KG_dir]}/super_pop_keep.list "
        ") &> {log}"


rule all_sblup_prep:
    input: 
        expand("{}/sblup/{{study.study_id}}/1KGPhase3.w_hm3.{{study.ancestry}}.{{study.study_id}}".format(config['Base_sumstats_dir']), study=studies.itertuples())
