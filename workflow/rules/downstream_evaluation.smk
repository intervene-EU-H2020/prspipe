
# things to run after the main pipeline
# these things go beyond what is implemented in GenoPred
# the outputs of these scripts can serve as inputs for cross-biobank meta-analyses


# the rules in here will only run for existing outputs of the main pipeline
# the main pipeline should be run with the rule all_get_best_models_ext first!

from glob import glob
from snakemake.io import glob_wildcards
import itertools

####################
# helper functions #
####################


def get_gz_pheno(wildcards):
    # handle gzipped input files
    infile = f'custom_input/{wildcards.bbid}/phenotypes/{wildcards.phenotype}.tsv.gz'
    if os.path.isfile(infile):
        return infile
    elif os.path.isfile(infile[:-3]):
        return infile[:-3]
    else:
        raise ValueError(f"Can't find phenotype file for phenotype '{wildcards.phenotype}'")


def file_not_empty(x):
    with open(x, 'r') as infile:
        for i, _ in enumerate(infile):
            if i > 1:
                return True # file has more than one row
        return False # file only has a header row


def join_comma(x):
    if isinstance(x,str):
        return x
    else:
        return ','.join(x)


def request_all(file_pattern, *args, **kwargs):

    phenotypes_config = [p.strip().split(',') for p in studies.name]
    study_to_phenotype_config = pd.DataFrame({'study_id':itertools.chain.from_iterable([s] * len(phenotypes_config[i]) for i, s in enumerate(studies.study_id)), 'phenotype':itertools.chain.from_iterable(phenotypes_config)})

    pe_files = []
    for i, row in study_to_phenotype_config.iterrows():
        for bbid in target_list.name.values:
            possible_pe_files = f'results/{bbid}/PRS_evaluation/{row.study_id}/{{superpop}}/{row.study_id}.{row.phenotype}.{{superpop}}.AllMethodComp.pred_eval.txt'
            possible_pe_files = expand(possible_pe_files, superpop=config['1kg_superpop'])
            for f in possible_pe_files:
                if not os.path.isfile(f):
                    continue
                else:
                    if file_not_empty(f):
                        pe_files.append(f)
    wc = glob_wildcards('results/{bbid}/PRS_evaluation/{study_ignore}/{superpop_ignore}/{study}.{pheno}.{superpop}.AllMethodComp.pred_eval.txt', pe_files)

    bsps = pd.DataFrame({'bbid':wc.bbid, 'study_id':wc.study, 'phenotype':wc.pheno, 'superpop':wc.superpop})

    files = expand(file_pattern, zip, bbid=bsps.bbid, study=bsps.study_id, phenotype=bsps.phenotype, superpop=bsps.superpop, allow_missing=True)

    if kwargs:
        files = expand(files, **kwargs)

    return files



##############################
# Calculate biobank MultiPRS #
##############################

rule predict_multiprs:
    # predict MultiPRS using weights and pre-computed standardized scores (profiles)
    input:
        model_rds='results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{phenotype}.{superpop}.AllMethodComp.{group}.model.rds'
    output:
        weights='results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{phenotype}.{superpop}.AllMethodComp.{group}.multiprs_weights.tsv',
        pgs='results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{phenotype}.{superpop}.AllMethodComp.{group}.multiprs_profiles.tsv.gz'
    log:
        'logs/predict_multiprs/{study}/{study}.{bbid}.{phenotype}.{superpop}.{group}.log'
    threads:
        1
    resources:
        mem_mb=32000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home",
        time="03:00:00",
        partition='cpu'
    singularity:
        config['singularity']['all']
    shell:
        "("
        "Rscript workflow/scripts/R/predict_multiprs.R "
        "--model_rds {input[model_rds]} "
        ") &> {log} "


eval_methods = list(M for M in config['prs_methods'] + ['All'] if M not in ['dbslmm','sbayesr']) # list of methods with more than one score
groups = [m.replace('_','.') +'_group' for m in eval_methods]

rule all_predict_multiprs:
    input:
        request_all(rules.predict_multiprs.output.weights, group=groups),
        request_all(rules.predict_multiprs.output.pgs, group=groups)


rule merge_multiprs:
    input:
        pgs = request_all(rules.predict_multiprs.output.pgs, group=groups),
        weights = request_all(rules.predict_multiprs.output.weights, group=groups)
    output:
        merged_profiles = request_all('results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{phenotype}.{superpop}.MultiPRS_profiles.merged.txt.gz'),
        merged_weights = request_all('results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{phenotype}.{superpop}.MultiPRS_weights.merged.tsv.gz')
    threads:
        1
    resources:
        mem_mb=64000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home",
        time="03:00:00",
        partition='cpu'
    singularity:
        config['singularity']['all']
    script:
        'workflow/scripts/R/merge_multiprs.R'



##############################################################
# calculate score-score correlations and performance metrics #
##############################################################


rule metrics_and_scor_train_test:
    # exports metrics and score-score correlations split by training/test
    # expect 2-3h runtime for 450k samples and 480 scores
    input:
        profiles = expand('results/{bbid}/prs/{method}/{study}/{superpop}/1KGPhase3.w_hm3.{study}.{superpop}.profiles', method = config['prs_methods'], allow_missing=True),
        biobank_multiprs = 'results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{phenotype}.{superpop}.MultiPRS_profiles.merged.txt.gz',
        pheno = get_gz_pheno,
        train_ind = 'results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{phenotype}.{superpop}.AllMethodComp.train_ind.txt.gz',
        test_ind = 'results/{bbid}/PRS_evaluation/{study}/{superpop}/{study}.{phenotype}.{superpop}.AllMethodComp.test_ind.txt.gz'
    output:
        mean_sd = 'results/{bbid}/metrics_and_scor_train_test/{study}/{study}.{phenotype}.{superpop}.mean_sd.tsv',
        metrics = 'results/{bbid}/metrics_and_scor_train_test/{study}/{study}.{phenotype}.{superpop}.metrics.tsv',
        scor = 'results/{bbid}/metrics_and_scor_train_test/{study}/{study}.{phenotype}.{superpop}.scor.tsv'
    params:
        profiles = lambda wc, input: join_comma(input['profiles']),
        train_test = lambda wc, input: input['train_ind'].replace('.train_ind.txt.gz',''),
        out_prefix = lambda wc, output: output['scor'].replace('.scor.tsv','')
    log:
        'logs/metrics_and_scor_train_test/{bbid}/{study}/{bbid}.{study}.{phenotype}.{superpop}.log'
    threads:
        1
    resources:
        mem_mb=128000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home",
        time="16:00:00",
        partition='cpu'
    singularity:
        config['singularity']['all']
    shell:
        "("
        "Rscript workflow/scripts/R/metrics_and_score_correl.R "
        "--profiles {params[profiles]} "
        "--train_test {params[train_test]} "
        "--pheno {input[pheno]} "
        "--out {params[out_prefix]} "
        "--boot 1000 "
        "--sample_cor 50000 "
        "--debug "
        ") &> {log} "


rule all_metrics_and_scor_train_test:
    input:
        request_all('results/{bbid}/metrics_and_scor_train_test/{study}/{study}.{phenotype}.{superpop}.mean_sd.tsv'),
        request_all('results/{bbid}/metrics_and_scor_train_test/{study}/{study}.{phenotype}.{superpop}.metrics.tsv'),
        request_all('results/{bbid}/metrics_and_scor_train_test/{study}/{study}.{phenotype}.{superpop}.scor.tsv')


rule metrics_and_scor_full:
    # exports metrics and score-score correlations for full sample
    # expect 2-3h runtime for 450k samples and 480 scores
    input:
        profiles = ancient(expand('results/{bbid}/prs/{method}/{study}/{superpop}/1KGPhase3.w_hm3.{study}.{superpop}.profiles',
                           method = config['prs_methods'], allow_missing=True)),
        pheno = ancient(get_gz_pheno)
    output:
        mean_sd = 'results/{bbid}/metrics_and_scor_full/{study}/{study}.{phenotype}.{superpop}.mean_sd.tsv',
        metrics = 'results/{bbid}/metrics_and_scor_full/{study}/{study}.{phenotype}.{superpop}.metrics.tsv'
    params:
        profiles = lambda wc, input: join_comma(input['profiles']),
        out_prefix = lambda wc, output: output['metrics'].replace('.metrics.tsv','')
    log:
        'logs/metrics_and_scor_full/{bbid}/{study}/{bbid}.{study}.{phenotype}.{superpop}.log'
    threads:
        1
    resources:
        mem_mb=128000,
        misc="--container-image=/dhc/groups/intervene/prspipe_0_1_1.sqsh --no-container-mount-home",
        time="16:00:00",
        partition='cpu'
    singularity:
        config['singularity']['all']
    shell:
        "("
        "Rscript workflow/scripts/R/metrics_and_score_correl.R "
        "--profiles {params[profiles]} "
        "--pheno {input[pheno]} "
        "--out {params[out_prefix]} "
        "--boot 1000 "
        "--calc_cor F "
        "--debug "
        ") &> {log} "


rule all_metrics_and_scor_full:
    input:
        request_all(rules.metrics_and_scor_full.output)
