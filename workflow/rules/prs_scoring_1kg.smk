
####################################
# Polygenic scoring of 1kg samples #
####################################


rule prs_scoring_1kg_ancestry:
    input:
        genotypes = expand('resources/1kg/1KGPhase3.w_hm3.chr{chrom}.{ext}', chrom=range(1,23), ext=['bed','bim','fam']),
        keep = 'resources/1kg/keep_files/{ancestry}_samples.keep',
        score = 'prs/{method}/{study}/1KGPhase3.w_hm3.{study}.score.gz',
        afreq = expand('resources/1kg/freq_files/{ancestry}/1KGPhase3.w_hm3.{ancestry}.chr{chrom}.afreq', chrom=range(1,23), allow_missing=True)
    output:
        full = 'prs/{method}/{study}/1KGPhase3.w_hm3.{ancestry}.{study}.profiles.tsv.gz',
        hla = 'prs/{method}/{study}/1KGPhase3.w_hm3.{ancestry}.{study}.hla_profiles.tsv.gz'
    log:
        'logs/prs_scoring_1kg_ancestry/{study}/{ancestry}/{method}.log'
    resources:
        threads=1,
        mem_mb=16000,
        partition='hpcpu',
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home"
    shell:
         "("
         "Rscript workflow/scripts/R/simple_per_chromosome_scorer.R "
         "--target_chr resources/1kg/1KGPhase3.w_hm3.chr "
         "--keep {input[keep]} "
         "--score {input[score]} "
         "--afreq resources/1kg/freq_files/{wildcards[ancestry]}/1KGPhase3.w_hm3.{wildcards[ancestry]}.chr "
         "--out prs/{wildcards[method]}/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[ancestry]}.{wildcards[study]} "
         "--mem {resources[mem_mb]} "
         "--plink2 bin/plink2 && "
         "Rscript workflow/scripts/R/simple_per_chromosome_scorer.R "
         "--target_chr resources/1kg/1KGPhase3.w_hm3.chr "
         "--keep {input[keep]} "
         "--score {input[score]} "
         "--afreq resources/1kg/freq_files/{wildcards[ancestry]}/1KGPhase3.w_hm3.{wildcards[ancestry]}.chr "
         "--out prs/{wildcards[method]}/{wildcards[study]}/1KGPhase3.w_hm3.{wildcards[ancestry]}.{wildcards[study]} "
         "--plink2 bin/plink2 "
         "--mem {resources[mem_mb]} "
         "--hla_only TRUE "
         ") &> {log}"


rule all_prs_scoring_1kg_ancestry_superpop:
    input:
        expand(rules.prs_scoring_1kg_ancestry.output, method=config['prs_methods'], study=studies['study_id'], ancestry=config['1kg_superpop'])
        
        
localrules:
    all_prs_scoring_1kg_ancestry_superpop