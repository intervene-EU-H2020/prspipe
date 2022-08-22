
# https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html
 
import itertools

bbids = target_list['name'].values

############################
# Harmonization to HapMap3 #
############################

# See rules in genotype_harmonizatin.smk (harmonize_target_genotypes)

####################
# Ancestry scoring #
####################

# needs: 1000G genotypes, keep files, and integrated_call_samples_v3.20130502.ALL.panel_small

rule ancestry_scoring_ext:
    # ancestry scoring of custom genotypes
    # recommended 15G of RAM
    input:
        super_pop_keep=rules.create_ancestry.output['super_pop_keep'],
        harmonised_geno=expand("custom_input/{bbid}/genotypes/chr{chr}.bed", chr=range(1,23), allow_missing=True)
    output:
        eigenvec_sp=expand('results/{{bbid}}/Ancestry_identifier/AllAncestry.{superpop}.eigenvec', superpop=config['1kg_superpop']),
        keep_sp=expand('results/{{bbid}}/Ancestry_identifier/AllAncestry.{superpop}.keep', superpop=config['1kg_superpop']),
        scale_sp=expand('results/{{bbid}}/Ancestry_identifier/AllAncestry.{superpop}.scale', superpop=config['1kg_superpop']),
        eigenvec='results/{bbid}/Ancestry_identifier/AllAncestry.eigenvec',
        eigenvec_var='results/{bbid}/Ancestry_identifier/AllAncestry.eigenvec.var',
        log='results/{bbid}/Ancestry_identifier/AllAncestry.log',
        model_pred_keep=expand('results/{{bbid}}/Ancestry_identifier/AllAncestry.model_pred.{superpop}.keep', superpop=config['1kg_superpop']),
        model_pred='results/{bbid}/Ancestry_identifier/AllAncestry.model_pred',
        pop_model='results/{bbid}/Ancestry_identifier/AllAncestry.pop_model.rds',
        model_details='results/{bbid}/Ancestry_identifier/AllAncestry.pop_model_prediction_details.txt',
        scale='results/{bbid}/Ancestry_identifier/AllAncestry.scale'
    log:
        "logs/ancestry_scoring_ext/{bbid}.log"
    singularity:
        config['singularity']['all']
    resources:
        mem_mb=64000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh --no-container-mount-home",
        time="12:00:00"
    shell:
        "( "
        "{config[Rscript]} {config[GenoPred_dir]}/Scripts/Ancestry_identifier/Ancestry_identifier.R "
        "--target_plink_chr custom_input/{wildcards[bbid]}/genotypes/chr "
        "--ref_plink_chr resources/1kg/1KGPhase3.w_hm3.chr "
        "--n_pcs 100 "
        "--plink {config[plink1_9]} "
        "--plink2 {config[plink2]} "
        "--output results/{wildcards[bbid]}/Ancestry_identifier/AllAncestry "
        "--ref_pop_scale {input[super_pop_keep]} "
        "--pop_data resources/1kg/integrated_call_samples_v3.20130502.ALL.panel_small "
        ") &> {log} "

    
rule ancestry_outlier_ext:
    input:
        rules.ancestry_scoring_ext.output,
        expand(rules.harmonize_target_genotypes.output, chr=range(1,23), allow_missing=True),
        keep_files = rules.ancestry_scoring_ext.output.model_pred_keep
    output:
        # many other output files not listed here.
        touch('results/{bbid}/Ancestry_identifier/outlier_detection/ok'),
        keep_list = "results/{bbid}/Ancestry_identifier/outlier_detection/AllAncestry.model_pred.keep_list",
        keep_files = expand("results/{{bbid}}/Ancestry_identifier/outlier_detection/AllAncestry.QC.{superpop}.keep", superpop=config['1kg_superpop'])
    resources:
        mem_mb=20000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh --no-container-mount-home",
        time="12:00:00"
    params:
        keep_files = lambda wc, input: ','.join(input['keep_files'])
    log:
        'logs/ancestry_outlier_ext/{bbid}.log'
    shell:
        # "ls results/{wildcards[bbid]}/Ancestry_identifier/*model_pred.*.keep > {output[keep_list]} && "
        "( "
        "echo {params[keep_files]} | tr ',' '\\n' > {output[keep_list]} && "
        "{config[Rscript]} {config[GenoPred_dir]}/Scripts/Population_outlier/Population_outlier.R "
        "--target_plink_chr custom_input/{wildcards[bbid]}/genotypes/chr "
        "--target_keep {output[keep_list]} "
        "--n_pcs 8 "
        "--maf 0.05 "
        "--geno 0.02 "
        "--hwe 1e-6 "
        "--memory 15000 "
        "--plink {config[plink1_9]} "
        "--plink2 {config[plink2]} "
        "--output results/{wildcards[bbid]}/Ancestry_identifier/outlier_detection/AllAncestry.QC "
        ") &> {log} "
    
if config['enable_outlier_detection']:
    keep_file_pattern = "results/{bbid}/Ancestry_identifier/outlier_detection/AllAncestry.QC.{superpop}.keep"
else:
    keep_file_pattern = "results/{bbid}/Ancestry_identifier/AllAncestry.{superpop}.keep"

rule all_ancestry_scoring_ext:
    input:
        expand(keep_file_pattern, bbid=target_list.name.values, superpop=config['1kg_superpop'])

rule calculate_maf_ancestry_ext:
    # calculate allele frequencies for different superpopulations
    # NA-values are dropped (TODO: make sure this doesn't lead to bugs downstream)
    input:
        keep_sp=keep_file_pattern,
        harmonised_geno="custom_input/genotypes/{bbid}/chr{chr}.bed"
    output:
        afreq='results/{bbid}/Ancestry_identifier/AllAncestry.{superpop}.chr{chr}.frq'
    params:
        geno_prefix=lambda wc, input: input['harmonised_geno'][:-4],
        out_prefix=lambda wc, output: output['afreq'][:-4]
    log:
        'logs/calculate_maf_ancestry_ext/{bbid}/{superpop}_{chr}.log'
    singularity:
        config['singularity']['all']
    shell:
        '('
        '{config[plink1_9]} '
        '--bfile {params[geno_prefix]} '
        '--freq '
        '--keep {input[keep_sp]} '
        '--out {params[out_prefix]}_tmp ; '
        'awk \'$5 != "NA"\' {params[out_prefix]}_tmp.frq > {params[out_prefix]}.frq && rm {params[out_prefix]}_tmp.frq '
        ') &> {log} '


rule all_calculate_maf_ancestry_ext:
    # runs rule above for all ancestries and chromosomes
    # Note: this fails if some ancestries are not represented, which is likely. Use "run.sh --keep-going ..." to force execution for all represented ancestries
    input:
        expand(rules.calculate_maf_ancestry_ext.output, chr=range(1,23), superpop=config['1kg_superpop'], bbid=bbids)
        

#################
# 10K reference #
#################

# we skip the generation of a 10K reference as done in the original GenoPred paper, see (deprecated) ukbb.smk if we wished to implement this for external biobanks as well.

        
#############################>>
# START: Polygenic scoring  #>>
#############################>>

rule validate_setup_ext:
    # requests all necessary outputs for the rules below.
    input:
        expand('prs/{method}/{study}/ok', method=config['prs_methods'], study=studies.study_id),
        expand('prs/{method}/{study}/1KGPhase3.w_hm3.{study}.score.gz', method=config['prs_methods'], study=studies.study_id),
        expand('prs/{method}/{study}/1KGPhase3.w_hm3.{study}.{superpop}.scale', method=config['prs_methods'], study=studies.study_id, superpop=config['1kg_superpop'])


wildcard_constraints:
    score_id="[A-Za-z\\.0-9_\\-]+"


def get_mem_mb_scaled_polygenic_scorer(wildcards, input, attempt):
    # this thing actually doesn't appear to need much memory at all (?)
    # does more memory make it faster (?) ... I doubt it.
    # requests roughly 16G for 500,000 people
    i = 0
    with open(input['target_keep'], 'r') as infile:
        for j, _ in enumerate(infile):
            i += 1
    if i == 0:
        print('Warning, empty keep file: ' + input['keep_file'])
    return max(4000, int(i * 0.03 * 1.5**(attempt-1)))
    
    

rule run_scaled_polygenic_scorer:
    # general purpose rule to run scaled_polygenic_scorer
    input:
        score='prs/{method}/{study}/{score_id}.score.gz',
        scale='prs/{method}/{study}/{score_id}.{superpop}.scale',
        ref_freq_chr=expand('resources/1kg/freq_files/{{superpop}}/1KGPhase3.w_hm3.{{superpop}}.chr{chr}.afreq', chr=range(1,23)),
        target_plink=expand('custom_input/{{bbid}}/genotypes/chr{chr}.{suffix}', chr=range(1,23), suffix=['bed','bim','fam']),
        target_keep=keep_file_pattern
    output:
        ok=touch('results/{bbid}/prs/{method}/{study}/{superpop}/{score_id}.{superpop}.done'),
        profiles='results/{bbid}/prs/{method}/{study}/{superpop}/{score_id}.{superpop}.profiles'
    params:
        geno_prefix=lambda wc, input: input['target_plink'][0][:-5],
        freq_prefix=lambda wc, input: input['ref_freq_chr'][0][:-7],
        out_prefix=lambda wc, output: output['ok'].replace('.done','')
    resources:
        mem_mb=get_mem_mb_scaled_polygenic_scorer,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh --no-container-mount-home",
        time="03:00:00",
        partition='vcpu,hpcpu'
    log:
        "logs/run_scaled_polygenic_scorer/{bbid}/{study}/{method}_{score_id}.{superpop}.log"
    shell:
        "({config[Rscript]} {config[GenoPred_dir]}/Scripts/Scaled_polygenic_scorer/Scaled_polygenic_scorer_plink2.R "
        "--target_plink_chr {params.geno_prefix} "
        "--target_keep {input.target_keep} "
        "--ref_score {input.score} "
        "--ref_scale {input.scale} "
        "--ref_freq_chr {params.freq_prefix} "
        "--plink2 {config[plink2]} "
        "--pheno_name {wildcards.study} "
        "--output {params[out_prefix]} "
        "--freq_format plink2 "
        "--safe TRUE) &> {log}"


# run target scoring for all methods defined in the pipeline
rule all_run_scaled_polygenic_scorer:
    input:
        expand('results/{bbid}/prs/{method}/{study}/{superpop}/1KGPhase3.w_hm3.{study}.{superpop}.done',
               bbid=target_list.name,
               study=studies.study_id.values,
               superpop=config['1kg_superpop'],
               method=config['prs_methods'])


# run target scoring for pT + clump
# TODO: all super-populations
rule all_target_prs_pt_clump:
    input:
        expand('results/{bbid}/prs/{method}/{study}/{superpop}/1KGPhase3.w_hm3.{study}.{superpop}.done',
               bbid=target_list.name,
               study=studies.study_id.values,
               superpop=config['1kg_superpop'],
               method=['pt_clump'])
               

###########################################################
# target scoring for all files that match general pattern #
###########################################################

from glob import glob
import os

def collect_prs_outputs():

    # function that scans the prs/ directory for scores, and returns the corresponding output files to request in order to run predictions

    score_dirs = glob('prs/*/*/')
    bbids = target_list.name.values
    
    outfile_request = []

    for basedir in score_dirs:
        _, method, study, _ = basedir.split('/') 
    
        score_files = glob(basedir + '*score.gz')
        
        if not len(score_files):
            continue
    
        score_ids = [ x.split('/')[-1].replace('.score.gz','') for x in score_files ]
        scale_files = list( s for s in glob(basedir + '*.scale') if s.split('/')[-1].split('.')[-2] in ['AFR','AMR','EAS','EUR','SAS'] and '.'.join(s.split('/')[-1].split('.')[:-2]) in score_ids)
        
        if not len(scale_files):
            print('Warning: found score-files but no scale-files in directory {}'.format(basedir))
            continue
            
        for s in scale_files:
            
            superpop = s.split('/')[-1].split('.')[-2]
            score_id = '.'.join(s.split('/')[-1].split('.')[:-2])
            
            request = ['results/{bbid}/prs/{method}/{study}/{superpop}/{score_id}.{superpop}.profiles'.format(bbid=b, method=method, study=study, superpop=superpop, score_id=score_id) for b in bbids ]
            
            outfile_request += request
        
    return outfile_request
    

rule all_target_prs_available:
    input:
        collect_prs_outputs()
        
localrules: all_target_prs_available

############################>>
# START: Score evaluation  #>>
############################>>
    
rule model_eval_ext_prep:
    # preparation for model evaluation, requests for all studies and all methods specified in the config.yaml
    input:
        expand('results/{bbid}/prs/{method}/{study}/{superpop}/1KGPhase3.w_hm3.{study}.{superpop}.profiles', method=config['prs_methods'], allow_missing=True)
    output:
        predictors='results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.AllMethodComp.predictor_groups'
    run:
        # prepare a file with predictors, grouped by method
        with open(output['predictors'], 'w') as outfile:
            outfile.write('predictors group\n')
            for i in input:
                method = i.split('/')[3]
                outfile.write('{} {}\n'.format(i, method))
                
localrules: model_eval_ext_prep

def check_gz_pheno_tsv(wc):
    # checks for phenotype files and prints a warning if they are missing
    pheno = wc['pheno']
    infile = 'custom_input/{bbid}/phenotypes/{pheno}.tsv.gz'.format(bbid = wc.bbid, pheno = pheno)
    
    if os.path.isfile(infile):
        return infile
    elif os.path.isfile(infile[:-3]):
        return infile[:-3]
    else:
        print('Error: Could not find phenotype file for target data "{}" and phenotype "{}", should be either "{}" or "{}" '.format(wc.study, pheno, infile, infile[:-3]))
        return infile


def get_mem_mb_model_eval_ext(wildcards, input, attempt):
    # calculate peak memory allocated to job 
    i = 0
    with open(input['keep_file'], 'r') as infile:
        for j, _ in enumerate(infile):
            i += 1
    if i == 0:
        print('Warning, empty keep file: ' + input['keep_file'])
    return max(8000, int(i * 0.23 * 1.5**(attempt-1)))



rule model_eval_ext:
    # evaluate score performance using target phenotype data
    input:
        predictors = rules.model_eval_ext_prep.output.predictors,
        pheno_file = check_gz_pheno_tsv,
        keep_file = keep_file_pattern,
        profiles = rules.model_eval_ext_prep.input
    output:
        assoc='results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{pheno}.{superpop}.AllMethodComp.assoc.txt',
        pred_eval='results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{pheno}.{superpop}.AllMethodComp.pred_eval.txt'
        #train_ind='results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{pheno}.{superpop}.AllMethodComp.train_ind.txt.gz',
        #test_ind='results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{pheno}.{superpop}.AllMethodComp.test_ind.txt.gz'
    params:
        prev = lambda wc: prevalence[wc.pheno],
        out_prefix = lambda wc, output: output['assoc'].replace('.assoc.txt','')
    threads:
        4
    resources:
        mem_mb=get_mem_mb_model_eval_ext,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh --no-container-mount-home",
        time="16:00:00",
        partition='hpcpu,vcpu'
    log:
        'logs/model_eval_ext/{bbid}/{study}.{pheno}.{superpop}.log'
    singularity:
        config['singularity']['all']
    shell:
        "("
        "export OPENBLAS_NUM_THREADS=1; "
        "{config[Rscript]} {config[GenoPred_dir]}/Scripts/Model_builder/Model_builder_V2.R "
        "--pheno {input[pheno_file]} "
        "--predictors {input[predictors]} "
        "--n_core {threads} "
        "--compare_predictors T "
        "--eval_only T "
        "--assoc T "
        "--min_case_val 10 "
        "--outcome_pop_prev {params[prev]} "
        "--out {params[out_prefix]} "
        "--save_group_model T "
        "--keep {input[keep_file]} "
        ") &> {log}"


def get_available_phenotypes(bbid):
    # return all available study (i.e., GWAS) to phenotype combinations for a specific biobank (bbid)
    pheno_path = glob('custom_input/{bbid}/phenotypes/*.tsv.gz'.format(bbid = bbid))
    pheno_path += glob('custom_input/{bbid}/phenotypes/*.tsv'.format(bbid = bbid))

    if not len(pheno_path):
        print('Warning, no phenotypes found for target data "{}"'.format(bbid))
        return [], []
    
    phenos = [ p.split('/')[-1] for p in pheno_path ]
    phenos = [ p[:-4] if p.endswith('tsv') else p[:-7] for p in phenos ]


    study_list = []
    pheno_list = []

    for _, row in studies.iterrows():
        study_phenos = row['name'].split(',')
        
        for p in study_phenos:
            if p in phenos:
                study_list += [row['study_id']]
                pheno_list += [p]

    return study_list, pheno_list


def get_possible_eval_outputs_ancestry(wc):
    # input function for a rule below
    # gets files for a single superpopulation
    
    study_list, pheno_list  = get_available_phenotypes(wc.bbid)

    superpop = [wc.superpop] * len(study_list)
    bbid = [bbid] * len(study_list)

    return expand(rules.model_eval_ext.output, zip, bbid=bbid, superpop=superpop, pheno=pheno_list, study = study_list)


def get_possible_eval_outputs(wc):
    # input function for a rule below
    # gets files for all superpopulations
    
    study_list, pheno_list  = get_available_phenotypes(wc.bbid)
    n_eval = len(study_list)

    study_list = study_list * len(config['1kg_superpop'])
    pheno_list = pheno_list * len(config['1kg_superpop'])

    superpop = []
    for s in config['1kg_superpop']:
        superpop += [s] * n_eval

    bbid = [wc.bbid] * len(study_list)

    return expand(rules.model_eval_ext.output, zip, bbid=bbid, superpop=superpop, pheno=pheno_list, study = study_list)


def get_possible_eval_ouputs_pred_eval(wc):
    outputs = get_possible_eval_outputs(wc)
    return list((o for o in outputs if '.pred_eval.' in o))


rule all_ancestry_model_eval_ext:
    # run score evaluation for all methods and all phenotypes for a specific superpopulation and target data
    input:
        get_possible_eval_outputs_ancestry
    output:
        touch('results/{bbid}/PRS_evaluation/{superpop}.ok')


rule all_target_model_eval_ext:
    # hack to make the input function work part 1
    input:
        get_possible_eval_outputs
    output:
        temp('results/{bbid}/PRS_evaluation/all_ancestry.ok')
    shell:
        'touch {output}'

rule all_model_eval_ext:
    # hack to make the input function work part 2
    # run score evaluation for all phenotypes, superpopulations and target data
    input:
        expand(rules.all_target_model_eval_ext.output, bbid=bbids)



localrules:
    all_ancestry_model_eval_ext,
    all_target_model_eval_ext,
    all_model_eval_ext
        
        
rule get_best_models_ext:
    # exctract the best performing models and their performance into a neat TSV for a single phenotype / biobank
    # it will often make more sense to run the rules below instead.
    input:
        pred_eval=rules.model_eval_ext.output['pred_eval']
    output:
        best_models_tsv='results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{pheno}.{superpop}.AllMethodComp.best_models.tsv'
    log:
        'logs/get_best_models_ext/{bbid}/{study}.{pheno}_{superpop}.log'
    singularity:
        config['singularity']['all']
    resources:
        mem_mb=8000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh --no-container-mount-home",
        time="00:10:00"
    shell:
        "("
        "{config[Rscript]} workflow/scripts/R/get_best_models_from_pred_eval.R "
        "--pred_eval {input[pred_eval]} "
        "--drop TRUE "
        "--out_prefix '{wildcards[study]}.{wildcards[pheno]}.{wildcards[superpop]}.AllMethodComp.best_models' "
        ") &> {log}"
        

rule biobank_get_best_models_ext:
    # extract the best performing models and their performance for all phenotypes for a single biobank
    # makes more sense to run this way because it doesn't use much compute
    input:
        get_possible_eval_ouputs_pred_eval
    output:
        touch('results/{bbid}/PRS_evaluation/all_get_best_models.ok')
    log:
        'logs/biobank_get_best_models_ext/{bbid}.log'
    singularity:
        config['singularity']['all']
    params:
        infiles=lambda wc, input: ' '.join(input)
    resources:
        mem_mb=8000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh --no-container-mount-home",
        time="00:30:00"
    shell:
        '('
        'for infile in {params[infiles]}; do '
        'BASENAME="$(basename ${{infile}})"; '
        'PREFIX="${{BASENAME%%.pred_eval.txt}}"; '
        '{config[Rscript]} workflow/scripts/R/get_best_models_from_pred_eval.R '
        '--pred_eval ${{infile}} '
        '--out_prefix ${{PREFIX}}.best_models '
        '--drop TRUE; '
        'done '
        ') &> {log}'


rule all_get_best_models_ext:
    # extract the best performing models for all phenotypes and biobanks
    input:
        expand(rules.biobank_get_best_models_ext.output, bbid=bbids)


localrules:
    all_get_best_models_ext
    
    
##################################
# Export best scores for sharing #
##################################


def get_existing_group_files():
    
    # finds .predictor_group file (only one needed for all ancestries)
    
    files = []
    
    for bbid in bbids:
        for study in studies.study_id.values:
            for ancestry in config['1kg_superpop']:
                path = 'results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.AllMethodComp.predictor_groups'.format(bbid=bbid, study=study, superpop=ancestry)
                if os.path.isfile(path):
                    files += [path]
                    break
    
    return files
                    
    

rule all_export_best_scores:
    # export best scores
    # essentially make a slimmed down copy of the .score.gz and .scale files in prs/..., based on the evaluation metrics of choice.
    # exports them to temp/prs/... the tmp/prs directory can then be shared, e.g. as a tar archive.
    # runs in a loop instead of separate jobs. The script could also be used in single jobs though (see below).
    # TODO: make this dynamically check for complete input files too
    params:
        group_files=get_existing_group_files(),
        ancestries=','.join(config['1kg_superpop'])
    output:
        touch('temp/all_export_best_scores.ok')
    singularity:
        config['singularity']['all']
    log:
        'logs/all_export_best_scores.log'
    resources:
        mem_mb=8000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_0.sqsh --no-container-mount-home",
        time="1:00:00"
    shell:
        '('
        'for infile in {params[group_files]}; do '
        '{config[Rscript]} workflow/scripts/R/export_best_scores.R '
        '--groups $infile '
        '--metric CrossVal_R '
        '--decreasing TRUE '
        '--N_max 20 '
        '--phenotype all '
        '--ancestries {params[ancestries]}; '
        'done '
        ') &> {log}'



# rule model_eval_custom_ext:
#     # Run Model_builder_V2.R for custom group - phenotype combinations
#     input:
#         predictors = 'custom_input/{bbid}/predictor_groups/{group}.predictor_groups',
#         pheno_file = 'custom_input/{bbid}/phenotypes/{pheno}.txt'
#     output:
#         assoc='results/{bbid}/PRS_evaluation_custom/{group}/{pheno}/AllMethodComp.assoc.txt',
#         pred_comp='results/{bbid}/PRS_evaluation_custom/{group}/{pheno}/AllMethodComp.pred_comp.txt',
#         pred_eval='results/{bbid}/PRS_evaluation_custom/{group}/{pheno}/AllMethodComp.pred_eval.txt'
#     params:
#         prev = lambda wc: prevalence[wc['pheno']],
#         out_prefix = lambda wc, output: output['assoc'].replace('.assoc.txt','')
#     threads:
#         16
#     log:
#         'logs/model_eval_custom_ext/{bbid}/{group}_{pheno}.log'
#     singularity:
#         config['singularity']['all']
#     shell:
#         "("
#         "{config[Rscript]} {config[GenoPred_dir]}/Scripts/Model_builder/Model_builder_V2.R "
#         "--pheno {input[pheno_file]} "
#         "--predictors {input[predictors]} "
#         "--n_core {threads} "
#         "--compare_predictors T "
#         "--assoc T "
#         "--outcome_pop_prev {params[prev]} "
#         "--out {params[out_prefix]} "
#         "--save_group_model T "
#         ") &> {log}"


# rule get_best_models_custom_ext:
#     # exctract the best performing models and their performance into a neat TSV
#     input:
#         pred_eval=rules.model_eval_custom_ext.output['pred_eval']
#     output:
#         best_models_tsv='results/{bbid}/PRS_evaluation_custom/{group}/{pheno}/best_models.tsv'
#     log:
#         'logs/get_best_models_custom/{bbid}/{group}_{pheno}.log'
#     singularity:
#         config['singularity']['all']
#     shell:
#         "("
#         "{config[Rscript]} workflow/scripts/R/get_best_models_from_pred_eval.R "
#         "--pred_eval {input[pred_eval]} "
#         "--drop TRUE "
#         ") &> {log} "
#     
