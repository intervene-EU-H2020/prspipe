
rule ldsc_download_1kg_genetic_map:
    # Download genetic map (and other things)
    output:
        "resources/1kg/genetic_map/1000GP_Phase3.tgz"
    shell:
        "wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz -O {output}"


rule ldsc_unpack_1kg_genetic_map:
    # extract the genetic map
    input:
        "resources/1kg/genetic_map/1000GP_Phase3.tgz"
    output:
        expand("resources/1kg/genetic_map/1000GP_Phase3/genetic_map_chr{chrom}_combined_b37.txt", chrom=range(1,23))
    shell:
        "cd resources/1kg/genetic_map/ && "
        "if [ -d 1000GP_Phase3 ]; then rm -r 1000GP_Phase3; fi; "
        "tar -xvf 1000GP_Phase3.tgz "


rule ldsc_cm_map_1kg:
    # rename and annotate variants with genetic map
    input:
        expand('resources/1kg/1KGPhase3.chr{chrom}.{s}', chrom=range(1,23), s=['bed','bim','fam']),
        expand("resources/1kg/genetic_map/1000GP_Phase3/genetic_map_chr{chrom}_combined_b37.txt", chrom=range(1,23))
    output:
        expand("resources/1kg/for_ldsc/1KGPhase3.chr{chrom}.{s}", chrom=range(1,23), s=['bed','bim','fam'])
    params:
        tmp_prefix="resources/1kg/for_ldsc/1KGPhase3.tmp.chr",
        out_prefix="resources/1kg/for_ldsc/1KGPhase3.chr"
    resources:
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    log:
        "logs/ldsc_cm_map_1kg.log"
    shell:
        "("
        "for chrom in {{1..22}}; do "
        "plink --cm-map resources/1kg/genetic_map/1000GP_Phase3/genetic_map_chr${{chrom}}_combined_b37.txt ${{chrom}} --bfile resources/1kg/1KGPhase3.chr${{chrom}} --make-just-bim --out {params[tmp_prefix]}${{chrom}}; "
        "awk -F \"\\t\" 'BEGIN{{OFS=\"\\t\"}}{{if ($2 == \".\"){{ $2 = $1 \":\" $4 \":\" $5 \":\" $6}} else{{$2 = $2 \":\" $4 \":\" $5 \":\" $6}}; print($0)}}' {params[tmp_prefix]}${{chrom}}.bim > {params[out_prefix]}${{chrom}}.bim && rm {params[tmp_prefix]}${{chrom}}.bim; "
        "ln -s -r resources/1kg/1KGPhase3.chr${{chrom}}.bed {params[out_prefix]}${{chrom}}.bed; "
        "ln -s -r resources/1kg/1KGPhase3.chr${{chrom}}.fam {params[out_prefix]}${{chrom}}.fam; "
        "done; "
        ") &> {log}"


rule ldsc_filter_1kg_mincount2_ancestry_plink2:
    # get non-singleton/missing variants by ancestry
    # and remove duplicate variants
    input:
        rules.create_ancestry.output,
        rules.ldsc_cm_map_1kg.output
    output:
        variant_select=expand("resources/1kg/for_ldsc/{{ancestry}}/1KGPhase3.{{ancestry}}.chr{chr}.select", chr=range(1,23))
    wildcard_constraints:
        ancestry='[A-Z]+'
    params:
        tmp_prefix=lambda wildcards: "resources/1kg/for_ldsc/" + wildcards['ancestry'] + "/1KGPhase3." + wildcards['ancestry'] + ".tmp",
        out_prefix=lambda wildcards: "resources/1kg/for_ldsc/" + wildcards['ancestry'] + "/1KGPhase3." + wildcards['ancestry'] + ".chr"
    log:
        "logs/ldsc_filter_1kg_mincount2_ancestry_plink2/{ancestry}.log"
    resources:
        time="01:00:00",
        mem_mb=4000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    singularity:
        config['singularity']['all']
    shell:
        "("
        "for chrom in {{1..22}}; do "
        "{config[plink2]} --bfile resources/1kg/for_ldsc/1KGPhase3.chr${{chrom}} "
        "--keep resources/1kg/keep_files/{wildcards[ancestry]}_samples.keep "
        "--freq counts cols=\"reffreq,altfreq\" "
        "--rm-dup exclude-all "
        "--out {params[tmp_prefix]}; "
        "awk '{{if(NR > 1){{{{if($3 < $2){{if($3 > 1){{print($1)}}}}else{{if($2 > 1){{print($1)}}}}}}}}}}' {params[tmp_prefix]}.acount > {params[out_prefix]}${{chrom}}.select && rm {params[tmp_prefix]}.acount; "
        "done "
        ") &> {log} "
        
        
rule ldsc_create_temporary_ancestry_genotype_files:
    # create temporary bed/bim/fam-files to calculate LD-scores with (ancestry specific)
    input:
        select_file=expand("resources/1kg/for_ldsc/{{ancestry}}/1KGPhase3.{{ancestry}}.chr{chrom}.select", chrom=range(1,23)),
        keep_file="resources/1kg/keep_files/{ancestry}_samples.keep",
        genotypes=rules.ldsc_cm_map_1kg.output
    output:
        bed=expand("results/ldsc/{{ancestry}}/1KGPhase3.{{ancestry}}.chr{chrom}.bed", chrom=range(1,23)),
        bim=expand("results/ldsc/{{ancestry}}/1KGPhase3.{{ancestry}}.chr{chrom}.bim", chrom=range(1,23)),
        fam=expand("results/ldsc/{{ancestry}}/1KGPhase3.{{ancestry}}.chr{chrom}.fam", chrom=range(1,23))
    log:
        "logs/ldsc_create_temporary_ancestry_genotype_files/{ancestry}.log"
    shell:
        "("
        "for chrom in {{1..22}}; do "
        "{config[plink2]} --bfile resources/1kg/for_ldsc/1KGPhase3.chr${{chrom}} "
        "--extract resources/1kg/for_ldsc/{wildcards[ancestry]}/1KGPhase3.{wildcards[ancestry]}.chr${{chrom}}.select "
        "--keep {input[keep_file]} "
        "--make-bed "
        "--out results/ldsc/{wildcards[ancestry]}/1KGPhase3.{wildcards[ancestry]}.chr${{chrom}}; done "
        ") &> {log}"
        
        
rule ldsc_calculate_stratified_ld_scores_chr:
    input:
        bim='results/ldsc/{ancestry}/1KGPhase3.{ancestry}.chr{chrom}.bim',
        bed='results/ldsc/{ancestry}/1KGPhase3.{ancestry}.chr{chrom}.bed',
        fam='results/ldsc/{ancestry}/1KGPhase3.{ancestry}.chr{chrom}.fam',
        annot='custom_input/ldsc/{annot}.{chrom}.annot.gz'
    output:
        annot=temp('results/ldsc/{ancestry}/{annot}.{chrom}.annot.gz'),
        ldscore='results/ldsc/{ancestry}/{annot}.{chrom}.l2.ldscore.gz',
        M='results/ldsc/{ancestry}/{annot}.{chrom}.l2.M',
        M_5_50='results/ldsc/{ancestry}/{annot}.{chrom}.l2.M_5_50'
        # ldsc_log='results/ldsc/{ancestry}/{annot}.{chrom}.log'
    params:
        bfile=lambda wc, input: input['bed'][:-4],
        annot_prefix=lambda wc, output: output['annot'][:-3]
    log:
        "logs/ldsc_calculate_stratified_ld_scores/{ancestry}/{annot}.{chrom}.log"
    conda:
        "../envs/ldsc.yaml"
    resources:
        mem_mb=16000,
        time="4:30:00",
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    shell:
        "("
        "python workflow/scripts/Python/ldsc_bim_annot_intersect.py "
        "--annot {input[annot]} "
        "--bim {input[bim]} "
        "--out {params[annot_prefix]} && "
        "gzip {params[annot_prefix]} && "
        ""
        "python {config[LDSC_dir]}/ldsc.py "
        "--ld-wind-cm 1.0 "
        "--bfile {params[bfile]} "
        "--annot {output[annot]} "
        "--l2 "
        "--out results/ldsc/{wildcards[ancestry]}/{wildcards[annot]}.{wildcards[chrom]} "
        ") &> {log}"
        
        
        
rule ldsc_calculate_stratified_ld_scores_annot:
    input:
        expand(rules.ldsc_calculate_stratified_ld_scores_chr.output, chrom=range(1,23), allow_missing=True)
    output:
        touch('results/ldsc/{ancestry}/{annot}.ok')
        
        
localrules: ldsc_calculate_stratified_ld_scores_annot
        


