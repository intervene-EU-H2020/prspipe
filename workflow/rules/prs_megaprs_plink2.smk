
# Download LDAK
rule download_ldak:
  output:
    directory("workflow/scripts/ldak")
  shell:
    "mkdir -p workflow/scripts/ldak; "
    "wget --no-check-certificate -O workflow/scripts/ldak/ldak5.1.linux_.zip https://dougspeed.com/wp-content/uploads/ldak5.1.linux_.zip; "
    "unzip workflow/scripts/ldak/ldak5.1.linux_.zip -d workflow/scripts/ldak/; "
    "rm workflow/scripts/ldak/ldak5.1.linux_.zip"


# Download LDAK map data
rule download_ldak_map:
  output:
    directory("resources/ldak/ldak_map")
  shell:
    "mkdir -p resources/ldak/ldak_map; " 
    "wget --no-check-certificate -O resources/ldak/ldak_map/genetic_map_b37.zip https://www.dropbox.com/s/slchsd0uyd4hii8/genetic_map_b37.zip; "
    "unzip resources/ldak/ldak_map/genetic_map_b37.zip -d resources/ldak/ldak_map/; "
    "rm resources/ldak/ldak_map/genetic_map_b37.zip"


# Download LDAK bld snp annotations
rule download_ldak_bld:
  output:
    directory("resources/ldak/ldak_bld")
  shell:
    "mkdir -p resources/ldak/ldak_bld; "
    "wget --no-check-certificate -O resources/ldak/ldak_bld/bld.zip https://genetics.ghpc.au.dk/doug/bld.zip; "
    "unzip resources/ldak/ldak_bld/bld.zip -d resources/ldak/ldak_bld/; "
    "rm resources/ldak/ldak_bld/bld.zip"


# Download LDAK high ld regions file
rule download_ldak_highld:
  output:
    "resources/ldak/ldak_highld/highld.txt"
  shell:
    "mkdir -p resources/ldak/ldak_highld; "
    "wget --no-check-certificate -O resources/ldak/ldak_highld/highld.txt https://dougspeed.com/wp-content/uploads/highld.txt"


rule prs_scoring_megaprs:
    input:
        rules.download_ldak.output,
        rules.download_ldak_map.output,
        rules.download_ldak_bld.output,
        rules.download_ldak_highld.output,
        hm3_gw=rules.extract_hm3_gw.output,
        hm3_chr=expand(rules.extract_hm3.output, chr=range(1,23)),
        qc_stats=lambda wc: expand(rules.QC_sumstats.output, ancestry = studies.ancestry[studies.study_id == wc.study], allow_missing=True),
        super_pop_keep=rules.create_ancestry.output['super_pop_keep']
    output:
        touch('prs/megaprs/{study}/ok')
    params:
        study_ancestry=lambda wc: studies.ancestry[studies.study_id == wc['study']].iloc[0]
    resources: 
        mem_mb=20000,
        time="14:00:00",
        misc="--container-image=/dhc/groups/intervene/prspipe_0_0_1.sqsh --no-container-mount-home"
    log:
        "logs/prs_scoring_megaprs/{study}.log"
    threads:
        5
    shell:
        "("
        "Rscript {config[GenoPred_dir]}/Scripts/ldak_mega_prs/ldak_mega_prs.R "
        "--ref_plink resources/1kg/1KGPhase3.w_hm3.GW "
        "--ref_keep resources/1kg/keep_files/{params[study_ancestry]}_samples.keep "
        "--sumstats {input[qc_stats]} "
        "--plink1 {config[plink1_9]} "
        "--plink2 {config[plink2]} "
        "--ldak workflow/scripts/ldak/ldak5.1.linux "
        "--ldak_map resources/ldak/ldak_map/genetic_map_b37 "
        "--ldak_tag resources/ldak/ldak_bld "
        "--ldak_highld resources/ldak/ldak_highld/highld.txt "
        "--memory {resources[mem_mb]} "
        "--n_cores {threads} "
        "--output prs/megaprs/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[study]} "
        "--ref_pop_scale {input[super_pop_keep]} "
        ") &> {log}"

rule all_prs_scoring_megaprs:
    input:
        expand(rules.prs_scoring_megaprs.output, study=studies.study_id)
