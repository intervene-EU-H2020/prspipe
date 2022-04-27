# Snakemake workflow: prspipe

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/prspipe.svg?branch=master)](https://travis-ci.org/snakemake-workflows/prspipe)

Snakemake pipeline to run Polygenic Risk Score (PRS) prediction. Implements and extends the [GenoPred](https://github.com/opain/GenoPred) pipeline, i.e. a reference standardized framework for the prediction of PRS using different state of the art methods using summary statistics.

## Authors

* Remo Monti & Sophie Wharrie
* [GenoPred](https://github.com/opain/GenoPred): Oliver Pain 

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

# Tutorial for Collaborators

The full pipeline can roughly be devided into two stages: (1) Download of summary statistics, running PRS methods and predicting on the 1000 Genomes data (public stage), and (2) prediction and evaluation using the scores from step 1, which includes hyper-parameter tuning using cross validation. The second step uses private data e.g. from the UK Biobank or your biobank of interest. **We are currently evaluating running the second step in other biobanks**.

Follow the steps below to make the pipeline work with your data (i.e., the target data). For collaborators, we distribute pre-computed scoring files from the different PRS methods, which means most of the pipeline can be skipped. However, the rules to run the different PRS methods are available to everyone - feel free to try!

In the steps below, we will install all the dependencies to run polygenic scoring for your target data. We will then harmonize the target genotype data and perform ancestry scoring. Finally, we will run polygenic scoring and score evaluation. Evaluation requires the target phenotype data.

**Steps that need internet access are marked with :globe_with_meridians: and steps that require access to sensitive data are marked with :rotating_light:.**

## TL;DR

For the lazy / impatient / advanced users, here is a breakdown of what you have to do, including the most relevant commands

1. :globe_with_meridians: Clone this repository and switch to the base directory for all subsequent steps
2. :globe_with_meridians: Use the prspipe docker/singularity container to run snakemake (or install snakemake and R-packages manually)
3. :globe_with_meridians: run the setup script
    ```
    bash install_basics.sh
    ```
5. edit the software paths in `config/config.yaml`, if necessary
6. :globe_with_meridians: run the snakemake rules to download and pre-process the 1000 Genomes reference
    ```
    bash run.sh all_setup
    bash run.sh cleanup_after_setup
    ```
7.  :globe_with_meridians: Download pre-calculated scores for pruning & thresholding + clump, MegaPRS, PRScs for 10 phenotypes
    ```
    bash run.sh download_test_data
    ```
9.  Edit the sample-sheet for your target genetic data (`config/target_list.tsv`)
10. 🚨 Predict scores
    ```
    bash run.sh all_target_prs_available
    ```
The tutorial below covers the steps above (and more) in greater detail.

## Preface on Singularity and Docker :takeout_box:

Genopred, i.e., the repository this pipeline depends on, relies on R and dependency management with conda (i.e., python). The pipeline itself is run with [snakemake](https://snakemake.bitbucket.io). Snakemake comes with built-in support for Singularity containers. In theory, different steps of the pipeline (correspinding to different snakemake *rules*), can be run in different containers. However, this pipeline only relies on a single container [available on dockerhub](https://hub.docker.com/r/rmonti/prspipe). This container works both with Docker and Singularity.

```
# get the container with docker
docker pull rmonti/prspipe:0.1.0

# ... or with singularity
# this should work without root permissions
singularity pull docker://rmonti/prspipe:0.1.0
```
> ⚠️ To run population outlier detection, you need at least version 0.0.3

## Preface on Snakemake :snake:
 
 Snakemake will handle all the steps to get from a set of input files (typically defined in a [*samplesheet*](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/studies.tsv) or [*config-file*](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/config.yaml)) to a set of output files. The user requests specific output files, and snakemake will figure out how to produce them given previously defined [*rules*](https://github.com/intervene-EU-H2020/prspipe/blob/main/workflow/rules/). Snakemake has been build with HPC clusters in mind, which often rely on schedulers such as [slurm](https://en.wikipedia.org/wiki/Slurm_Workload_Manager). Snakemake will work with the scheduler to distribute the computational workload over many different hosts (in parallel, if possible). Running snakemake this way typically requires setting up HPC-cluster-specific configuration files (example shown in `slurm/config.yaml`).

However, snakemake can also be run interactively on a single host (i.e, a single server or virtual machine). This will be easier for most users to set up. To run the polygenic scoring of your target data, this setup will be sufficient. Therefore, **these instructions will handle the interactive case**.

### running snakemake rules

To run steps of the pipeline, the user can request specific output files on the command-line. For example

```
snakemake --cores 1 --snakefile workflow/Snakefile resources/HapMap3_snplist/w_hm3.snplist
```

will request `resources/HapMap3_snplist/w_hm3.snplist`.

If the output files do not contain any [*wildcards*](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards), the user can also request to run *rules* using the rule name directly. Rules are the different steps defining the snakemake workflow and are located in `workflow/Snakefile` and `workflow/rules/...smk`.

For example, the file above could also be requested by the rule name:

```
# "download_hapmap3_snplist" is the rule name
snakemake --cores 1 --snakefile workflow/Snakefile download_hapmap3_snplist
```

the corresponding rule looks like this (located in [`workflow/rules/setup.smk`](https://github.com/intervene-EU-H2020/prspipe/blob/main/workflow/rules/setup.smk):

```
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
```

The `shell:`-directive defines what will be executed. As you can see, this rule also defines a log-file. The log files often contain useful diagnostic messages. If a rule fails, make sure to check both the snakemake output and the log file!

In order to avoid writing a long command every time you run snakemake, the shell script [`run.sh`](https://github.com/intervene-EU-H2020/prspipe/blob/main/run.sh) in the main directory can be used to wrap default parameters, i.e.

```
bash run.sh download_hapmap3_snplist
```
is equivalent to the two commands in the examples above.

## :globe_with_meridians: Basic Setup and Dependencies

Clone the repository. The base directory of the repository will be your working directory.

```
git clone git@github.com:intervene-EU-H2020/prspipe.git && cd prspipe
```

Harmonizing target genotype data and performing polygenic scoring require a number of software dependencies, which can be installed by running the `install_basics.sh`-script.

```
bash install_basics.sh
```

This will clone our fork of the [GenoPred repository](https://github.com/intervene-EU-H2020/GenoPred), download `qctool2` and the latest versions of `plink1` and `plink2`. If you already have these tools installed, you can also change the paths in `config/config.yaml`. However, this can lead to bugs :cockroach: and is not recommended.  The default software paths in `config/config.yaml` should work when using the pipeline container. 

> :warning: If `workflow/scripts/GenoPred` does not exist after running `bash install_basics.sh`, try to clone the repo manually: `git clone git@github.com:intervene-EU-H2020/GenoPred.git workflow/scripts/GenoPred`. Getting an error? Make sure you have access, and your git is [configured to use ssh](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)!

## :globe_with_meridians: Set up Snakemake using Singularity

Snakemake is installed in the prspipe container. If singularity is available, you can create a local container image which will be able to run the pipeline (at least the polygenic scoring and genotype harmonization steps):

```
# this will create a local image called prspipe.sif
singularity pull prspipe.sif docker://rmonti/prspipe:0.1.0
```

After that, you can run snakemake interactively by first starting an [interactive shell inside the container](https://sylabs.io/guides/3.7/user-guide/cli/singularity_shell.html), and then running snakemake:

```
# on the host, in the working directory:
singularity shell -c prspipe.sif

# inside the container:
snakemake --help

# should be the container's R (located at /usr/local/bin/R):
which R

# to exit the container
exit
```
> :warning: Note: When running with singularity make sure to clear environment variables related to R such as `R_LIBS`, `R_LIBS_SITE` and `R_LIBS_USER`.

Beware, if the path displayed after `which R` is *not* /usr/local/bin/R, your environment variables, `.bashrc` or `.bash_profile` might be interfering with R in the container. Try changing the parameter `Rscript:` inside [`config/config.yaml`](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/config.yaml#L32) to force using the container's R installation (see the comments in that file).

> Note: If you do **not** have singularity you can also [install R-packages manually](#manual-installation-of-r-packages). You will then also have to [install snakemake](#installing-snakemake-with-conda), and adjust software paths in `config/config.yaml`.

## Container mounts :file_folder:

:warning: When running the container, you have to "mount" any directories that you want access to, unless they are subdirectories of the main working directory. For example, if you want to access genotype data at `/some/random/path/`, you will have to mount this directory (i.e., make it accessible) inside the container. You can do this by specifying mounts on the command-line when running singularity. The command above becomes `singularity shell -c -B /some/random/ prspipe.sif`. This would make `/some/random` and all its sub-directories available inside the container.

## :globe_with_meridians: Download and process the 1000 Genomes Reference

Assuming you have snakemake running, you can now download and pre-process the 1000 Genomes data. This data will serve as an external reference throughout the analysis and is used for ancestry scoring. First we perform a "quiet" (`-q`) "dry-run" (`-n`).

```
export SNAKEMAKE_CORES=1
bash run.sh -n -q all_setup
```
The command above will display a summary of what will be done. To then actually run these steps:

```
bash run.sh all_setup
```
> :warning: Note: depending on your download speed and how many cores you have available this can take several hours. It will also potentially require a lot of disk-space (~90G). To save disk space, you can use a small number of cores (i.e., fewer parallel processes) in snakemake.

Once you have successfully completed these steps, you can clear up space by running

```
bash run.sh cleanup_after_setup
```

# Testing PRS methods with synthetic data
In the sections below, we will cover how to run the pipeline on synthetic data. This will help you verify your setup is working, and will show the basics of how the pipeline can be used. **this part can be skipped if you only want to run the score evaluation steps**

## Test your setup by running pruning & thresholding

The pipeline ships with synthetic data (genotypes, phenotypes and summary statistics) generated by Sophie Wharrie, Vishnu Raj and Zhiyu Yang. By default, the pipeline is configured to work with this synthetic data. Run the following rule, to extract the synthetic data:

```
bash run.sh initialize_synthetic
```

You can now create polygenic scores using plink and summary statistics from a GWAS run on synthetic data:

```
# first, check what will be executed:
bash run.sh -n -q prs/pt_clump/synth01/ok

# second, run the rule
bash run.sh prs/pt_clump/synth01/ok
```

`config/studies.tsv` by default is the sample-sheet that contains the information on available summary statistics:

| study_id | ancestry | n_cases | n_controls | ftp_address | local_path | binary | name |
| --- | --- | --- | --- | --- | --- | --- | --- |
| The GWAS acccession | The ancestry abbreviation (EUR, EAS, AMR, SAS and AFR) | The number of case samples | The number of control samples | The GWAS catalog ftp address of the ```.h.tsv.gz``` harmonized summary statistics file, given in the form ```/pub/databases/gwas/summary_statistics.../harmonised/....h.tsv.gz```) | alternatively a local path to a "munged" summary statistics file | phenotype is binary (yes/no) | name (descriptive name, like "T2D", can be repeated) |

[`config/studies.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/studies.tsv) is configured to work with summary statistics generated from synthetic data. In order to create scores with a specific PRS method "`{method}`" and GWAS stummary statistics `"{study}"`, the user can request output files which follow the pattern: `prs/{method}/{study}/ok`. Available methods are `dbslmm`,`lassosum`,`ldpred2`,`megaprs`,`prscs`,`pt_clump`.

>:warning:Note: This is just an example. Collaborators don't have to run these methods. **Their outputs will be distributed.** See steps below.

## Predict polygenic scores for synthetic dataset 

Check out the polygenic scoring files created by the rule above:

```
ls -1 prs/pt_clump/synth01/
```

>Output (only showing those that are common between methods):
>```
>1KGPhase3.w_hm3.synth01.AFR.scale
>1KGPhase3.w_hm3.synth01.AMR.scale
>1KGPhase3.w_hm3.synth01.EAS.scale
>1KGPhase3.w_hm3.synth01.EUR.scale
>1KGPhase3.w_hm3.synth01.SAS.scale
>1KGPhase3.w_hm3.synth01.log
>1KGPhase3.w_hm3.synth01.score.gz
>ok
>```
>While `1KGPhase3.w_hm3.synth01.score.gz` contains the weights for the different SNPs, the scale-files contain the mean and standard deviation of the scores for the different ancestries, which are used for normalization later. 

To predict all *available* scores for all ancestries and target data, run:

```
# dryrun
# this should only trigger scoring for the scores we just created above

bash run.sh -n -q all_target_prs_available
```

```
# run the scoring
bash run.sh all_target_prs_available
```

## Evaluate polygenic score performance for the synthetic dataset

We now want to evaluate the performance of the polygenic scores we just created. Because we only have scores for pruning & thresholding + clump (`pt_clump`) available, we edit `config/config.yaml` to only consider that method:

```
# this will swap the default
#     prs_methods: ['dbslmm','lassosum','ldpred2','megaprs','prscs','pt_clump']
# with
#     prs_methods: ['pt_clump']

sed -i -e "s/^prs_methods: .\+/prs_methods: ['pt_clump']/g" config/config.yaml
```
The synthetic phenotype data for phenotype "synthetic01" and target data "synth" are located at `custom_input/synth/phenotypes/synthetic01.tsv.gz`. In general, the pipeline will look for phenotype data in the directories matching the pattern `custom_input/{target-name}/phenotypes/{phenotype-name}.tsv.gz` (more about that later).

We can then run polygenic score evaluation and model selection by running the rule `all_model_eval_ext`

```
bash run.sh all_model_eval_ext
```

Finally, we revert our changes to the `config/config.yaml`:

```
sed -i "s/prs_methods: \['pt_clump'\]/prs_methods: ['dbslmm','lassosum','ldpred2','megaprs','prscs','pt_clump']/g" config/config.yaml
```
> :warning:Note: if the above sed-commands don't work for you, you can of course just manually edit config/config.yaml. To revert all changes, do `git checkout -- config/config.yaml`.


# Ancestry and polygenic scoring for new target genotype/phenotye data
The steps below will guide you through the process of setting up the pipeline to work with new target genotype/phenotype data. If your research environment does not have access to the internet and is located at a different location, you will have to make sure to transfer the entire working directory tree into your new environment. Also, you have to successfully run at least `bash run.sh all_setup` and `bash run.sh cleanup_after_setup` before transferring to the protected research environment. You will also have to transfer any conda environments you are using (if any), and the singularity/docker container image.

## :rotating_light: Setting up target genetic data
Target genetic data are defined using a tab-separated-values file (tsv), by default this is [`config/target_list.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/target_list.tsv). This file has three columns:

| name | path | type |
| --- | --- | --- |
| the name of the target data (e.g., "ukbiobank") | the path/prefix of the imputed per-chromosome genotype data (e.g.,"some/path/prefix_chr") | the format (see below) |

Genotype data paths are expected to follow the pattern `<prefix>_chr<1-22>.<.bed/.bim/.fam/.bgen/.vcf.gz>` (mind the gap, i.e., the underscore). For bgen-formatted data, the sample-file should be called called `<prefix>.sample`. Both hg19 and hg38 are accepted. Note that the path will have to be accessible inside the container, see the section on [container mounts](#container-mounts-file_folder) above.

Supported data types are 
- **samp_imp_plink1**: plink 1 formated data (.bed/.bim/.fam)
- **samp_imp_bgen_ref_first**, **samp_imp_bgen_ref_last** or **samp_imp_bgen_ref_unknown**: [Bgen](https://www.well.ox.ac.uk/~gav/bgen_format/) formatted data (.bgen/.sample). UK Biobank imputed genpotypes are ref-first, others typically ref-last. Specifying only "samp_imp_bgen" will default to ref-last.
- **samp_imp_vcf**: VCF format

As mentioned above, genotype data are expected to be split by chromosome. For plink 1 formatted data, the following bash one-liner will split data by chromosome. Asumming your main file is called `main.bed/.bim/.fam`:

```
# replace "main" with the prefix of your data
# this will generate files tmp_chr{1-22}.{bed,bim,fam}
for i in {1..22}; do plink --bfile main --chr $i --make-bed --out tmp_chr${i}; done
```

### :rotating_light: Target genotype harmonization

Prspipe uses the same script as [GenoPredPipe](https://github.com/opain/GenoPred/tree/master/GenoPredPipe) to convert all target genetic data to plink1 hard-calls. Genotypes are harmonized based on variant positions and allele codes to match the HapMap3 referece intersected with the 1000 Genomes variants. Any existing variant identifiers wil be replaced by rsids.

To perform harmonization of the target genetic data, run

```
bash run.sh all_harmonize_target_genotypes
```

Harmonization can take several hours to finish depending on the size of your data. Once completed, the harmonized target data are located at `custom_input/{target-name}/genotypes/`, where `target-name` matches the entry in the name-column of the `target_list.tsv` (see above).

## :rotating_light: Ancestry scoring

Once the target genotype data are harmonized, we can perform ancestry scoring:

```
bash run.sh all_ancestry_scoring_ext
```

this will create output files at `results/{target-name}/Ancestry_idenitfier/`. 

## :rotating_light: Setting up target phenotype data

> ⚠️ This section is under construction 🚧

Extract the phenotypes of interest. A list of phenotypes for initial testing will be distributed. How to extract phenotypes will be highly specific to the biobank you are working with and may include self-reported endpoints or data from patient records. Phenotypes should be placed in *separate* files with three tab-separated columns (no header): The family ID, the individual ID, and the Phenotype value (see the plink "pheno" format [here](https://www.cog-genomics.org/plink/1.9/input#pheno)), for example

|  |  |  |
| --- | --- | --- |
| fid_1 | iid_1 | 0 |
| fid_2 | iid_2 | 1 |
| ... | ... | ... |

An example is provided at [`resources/synthetic_data/pheno250.tsv.gz`](https://github.com/intervene-EU-H2020/prspipe/blob/main/resources/synthetic_data/). A script to convert the INTERVENE phenotype format to this format will be provided soon (any volunteers?). 

The IDs should have matches in the plink `.fam`-files generated above, but it is **not** necessary to have phenotype values for all genotyped individuals. This   way only a sub-set of individuals can be analysed for each phenotype to save computational effort. Also, if you do not want to analyse a specific phenotype, you can delete the row corresponding to that phenotype from [`config/studies.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/studies.tsv).

Place the phenotype files in `custom_input/{target-name}/phenotypes/`. Name them `{phenotype}.tsv`, where {phenotype} should match the entries in the "name"-column of [`config/studies.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/studies.tsv), i.e `./custom_input/ukbb/phenotypes/T2D.tsv`, `./custom_input/ukbb/phenotypes/BMI.tsv` and so on. **You can also gzip these files to save space.**

The pipeline does not need to be configured to look for files in `custom_input/phenotypes/*/`. `config/studies.tsv` can be left unchanged if you wish to analyse all phenotypes. Otherwise delete the rows corresponding to the phenotypes you do not wish to analyse from [`config/studies.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/studies.tsv).

### :rotating_light: Run polygenic risk scoring and evaluation for target data

> ⚠️ This section is under construction 🚧

Assuming you have downloaded pre-adjusted summary statistics as described above (`bash run.sh download_test_data`), you can now perform hyper-parameter tuning (model selection) on your data. First, perform a dryrun to check that the right jobs will be executed:
    
```
bash run.sh --dryrun --keep-going all_get_best_models_ext
```

the output should look something like this

```
[...]
Job counts:
        count   jobs
        ...
```
    if you see significantly more jobs listed, your pipeline is misconfigured. Please contact me in this case (remo.monti@hpi.de).

This will run ancestry scoring, identify individuals with EUR ancestry, and perform predictions and hyper-parameter tuning for those individuals. You can now run these steps with `bash run.sh --keep-going all_get_best_models_ext`.

--------------------------------------
    
# Other Documentation

Other information related to running the pipeline and setting it up for new summary statistics.

## Installing Snakemake with conda

To install snakemake with [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

> Note: make sure to update conda, or use mamba to get the latest version of snakemake i.e. `6.2.1`

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Manual installation of R packages

The pipeline relies heavily on [R](https://www.r-project.org/). There are two options to handle R-package dependencies: (1, **recommended**) running the pipeline using singularity containers, or (2) running with a local conda and R installation and installing R-packages manually.

If manually installing R-packages, the following commands will install all dependencies:

```
R
install.packages(c('data.table','doMC','optparse','foreach','caret','ggplot2','cowplot','glmnet','MLmetrics','e1071','stringr','verification', 'RcppArmadillo', 'Matrix', 'fdrtool', 'psych', 'bigsnpr', 'bigreadr', 'runonce', 'NbClust', 'GGally'), dependencies=TRUE)
    
install.packages("bin/lassosum/lassosum_0.4.5.tar.gz", repos=NULL, type="source")
```

## Incorporating new GWAS summary statistics

The code implements rules for automatically downloading and pre-processing summary statistics from the [GWAS catalog](https://www.ebi.ac.uk/gwas/home).

To configure and run the code for downloading summary statistics:
1. Select the studies from the GWAS catalog that you want to use
2. For each study, enter the details in the config/studies.tsv file:

| study_id | ancestry | n_cases | n_controls | ftp_address | local_path | binary | name |
| --- | --- | --- | --- | --- | --- | --- | --- |
| The study acccession number given in the GWAS catalog | The ancestry abbreviation (we currently support EUR, EAS, AMR, SAS and AFR) | The number of case samples | The number of control samples | The ftp address of the ```.h.tsv.gz``` harmonized summary statistics file, given in the form ```/pub/databases/gwas/summary_statistics.../harmonised/....h.tsv.gz```) | alternatively a local path to a "munged" summary statistics file | phenotype is binary (yes/no) | name |
| e.g. GCST000998 | e.g. EUR | e.g. 22233 | e.g. 64762 | e.g. /pub/databases/gwas/summary_statistics/GCST000001-GCST001000/GCST000998/harmonised/21378990-GCST000998-EFO_0001645.h.tsv.gz | e.g. ./munged_ss.tsv.gz | e.g. yes | e.g. CAD |

If summary statistics are not available in the harmonized format, consider using [this script](https://github.com/intervene-EU-H2020/prspipe/blob/workflow/scripts/R/munge_sumstats.R) to convert them to munged format.

3. Run the ```all_QC``` rule in the ```base_sumstats.smk``` rule file. Snakemake will download the summary statistics and run scripts to format and QC the data 

```
./run.sh all_QC
```

## Disease prevalence estimates
Model evaluation requires an estimate of disease prevalence in the general population for dichotomous traits. These can be found in the [`pop_prevalence.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/pop_prevalence.tsv). Make sure to add a corresponding row in that file when adding new endpoints. 
