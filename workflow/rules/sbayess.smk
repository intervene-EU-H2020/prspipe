
rule sbayess:
    input:
        sumstats = "resources/sumstats/{study}.EUR.cleaned.gz"
    output:
        sumstats = 'manuscript/polygenicity/{study}/{study}.sumstats.ma',
        imputedPerSnpN = 'manuscript/polygenicity/{study}/{study}.sbayess.imputedPerSnpN',
        Par = 'manuscript/polygenicity/{study}/{study}.sbayess.mcmcsamples.Par',
        snpRes = 'manuscript/polygenicity/{study}/{study}.sbayess.snpRes',
        SnpEffects = 'manuscript/polygenicity/{study}/{study}.sbayess.mcmcsamples.SnpEffects',
        parRes = 'manuscript/polygenicity/{study}/{study}.sbayess.parRes'
    log:
        'manuscript/polygenicity/{study}/{study}.log'
    resources:
        threads = 1,
        mem_mb = 48000,
        time = "12:00:00",
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home",
        partition='cpu'
    shell:
        "("
        "Rscript manuscript/polygenicity/run_sbayes_s.R "
        "--gctb bin/gctb/gctb "
        "--sumstats {input[sumstats]} "
        "--ld_matrix_chr resources/LD_matrix/sbayesr/UKBB/precomputed/EUR "
        "--chain_length 2000 "
        "--num_chains 4 "
        "--burn_in 200 "
        "--out manuscript/polygenicity/{wildcards[study]}/{wildcards[study]} "
        ") &> {log} "
        
        
rule sbayess_excl_mhc:
    input:
        sumstats = "resources/sumstats/{study}.EUR.cleaned.gz"
    output:
        sumstats = 'manuscript/polygenicity/{study}/{study}.excl_mhc.sumstats.ma',
        imputedPerSnpN = 'manuscript/polygenicity/{study}/{study}.excl_mhc.sbayess.imputedPerSnpN',
        Par = 'manuscript/polygenicity/{study}/{study}.excl_mhc.sbayess.mcmcsamples.Par',
        snpRes = 'manuscript/polygenicity/{study}/{study}.excl_mhc.sbayess.snpRes',
        SnpEffects = 'manuscript/polygenicity/{study}/{study}.excl_mhc.sbayess.mcmcsamples.SnpEffects',
        parRes = 'manuscript/polygenicity/{study}/{study}.excl_mhc.sbayess.parRes'
    log:
        'manuscript/polygenicity/{study}/{study}.excl_mhc.log'
    resources:
        threads = 1,
        mem_mb = 48000,
        time = "12:00:00",
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home",
        partition='hpcpu,vcpu'
    shell:
        "("
        "Rscript manuscript/polygenicity/run_sbayes_s.R "
        "--gctb bin/gctb/gctb "
        "--sumstats {input[sumstats]} "
        "--ld_matrix_chr resources/LD_matrix/sbayesr/UKBB/precomputed/EUR "
        "--chain_length 2000 "
        "--num_chains 4 "
        "--burn_in 200 "
        "--out manuscript/polygenicity/{wildcards[study]}/{wildcards[study]}.excl_mhc "
        "--exclude_mhc TRUE "
        ") &> {log} "
        
        
rule all_sbayess:
    input:
        expand(rules.sbayess.output, study=studies.study_id),
        expand(rules.sbayess_excl_mhc.output, study=studies.study_id),

        


