# Snakemake workflow: prspipe

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/prspipe.svg?branch=master)](https://travis-ci.org/snakemake-workflows/prspipe)

Snakemake pipeline to run Polygenic Risk Score (PRS) prediction. Implements and extends the [GenoPred](https://github.com/opain/GenoPred) pipeline, i.e. a reference standardized framework for the prediction of PRS using different state of the art methods using summary statistics.

## Authors

* Remo Monti & Sophie Wharrie
* [GenoPred](https://github.com/opain/GenoPred): Oliver Pain 

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

# Quick Start for Collaborators

The full pipeline can roughly be devided into two stages: (1) Download of summary statistics, running PRS methods and predicting on the 1000 Genomes data (public stage), and (2) prediction and evaluation using the scores from step 1, which includes hyper-parameter tuning using cross validation. The second step uses private data e.g. from the UK Biobank or your biobank of interest. **We are currently evaluating running the second step in other biobanks**.

Follow the steps below to make the pipeline work with your data (i.e., the target data) [`workflow/rules/external_biobanks.smk`](https://github.com/intervene-EU-H2020/prspipe/blob/main/workflow/rules/external_biobanks.smk). For collaborators, we distribute pre-computed scoring files from the different PRS methods, which means most of the pipeline can be skipped. However, the rules to run the different PRS methods are available to everyone - feel free to try!

In the steps below, we will install all the dependencies to run polygenic scoring for your target data. We will then harmonize the target genotype data. Finally, we will run polygenic scoring and score evaluation. Evaluation requires the target phenotype data.

Steps that need internet access are marked with :globe_with_meridians: and steps that require access to sensitive data are marked with :rotating_light:.

## Preface on Singularity and Docker :takeout_box:

Genopred, i.e., the repository this pipeline depends on, heavily relies on R and dependency management with conda (i.e., python). The pipeline itself is run with [snakemake](https://snakemake.bitbucket.io). Snakemake comes with built-in support for Singularity containers. In theory, different steps of the pipeline (correspinding to different snakemake *rules*), can be run in different containers. However, this pipeline only relies on a single container [available on dockerhub](https://hub.docker.com/r/rmonti/prspipe). This container works both with Docker and Singularity.

```
# get the container with docker
docker pull rmonti/prspipe:0.0.1

# ... or with singularity
singularity pull docker://rmonti/prspipe:0.0.1
```

## Preface on Snakemake :snake:

Snakemake reproducibly manages large workflows. Basically, it will handle all the steps to get from a set of input files (typically defined in a *samplesheet* or *config-file*) to a set of outout files. The user requests a set of output files, and snakemake will figure out how to produce them given previously defined *rules*. Snakemake has been build with HPC clusters in mind, which often rely on schedulers such as [slurm](https://en.wikipedia.org/wiki/Slurm_Workload_Manager). Snakemake will work with the scheduler to distribute the computational workload over many different hosts (in parallel, if possible). Running snakemake this way typically requires setting up HPC-cluster-specific configuration files (example shown in `slurm/config.yaml`).

However, snakemake can also be run interactively on a single host (i.e, a single server or virtual machine). This will be easier for most users to set up. To run the polygenic scoring of your target data, this setup will be sufficient. Therefore, these instructions will handle the interactive case.

To run steps of the pipeline, the user can request specific output files on the command-line. For example

```
snakemake --cores 1 --snakefile workflow/Snakefile resources/HapMap3_snplist/w_hm3.snplist
```

will request `resources/HapMap3_snplist/w_hm3.snplist`. There are many parameters that can be configured (see `snakemake --help`).

If the output file does not contain any *wildcards* (i.e., it is always called the same, no matter how the pipeline is configured), the user can also request to run specific *rules* using their name directly. Rules are the different steps defining the snakemake workflow and are located in `workflow/Snakefile` and `workflow/rules/...smk`.

For example, the file above could also be requested by the rule name:

```
snakemake --cores 1 --snakefile workflow/Snakefile download_hapmap3_snplist
```

the corresponding rule looks like this (located in `workflow/rules/setup.smk`):

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

As you can see, this rule also defines a log-file. The log files often contain useful diagnostic messages. If a rule fails, make sure to check both the snakemake output and the log file!

In order to not have to type a long command with snakemake parameters every time you run snakemake, the shell script `run.sh` in the main workflow directory can be used to wrap default parameters, i.e.

```
bash run.sh workflow/Snakefile download_hapmap3_snplist
```
is equivalent to the two commands in the examples above. 

## :globe_with_meridians: Basic Setup

Clone the repository. The base directory of the repository will be your working directory.

```
git clone git@github.com:intervene-EU-H2020/prspipe.git && cd prspipe
```

Harmonizing target genotype data and performing polygenic scoring require a number of software dependencies, which can be installed by running the `install_basics.sh`-script.

```
bash install_basics.sh
```

This will clone our fork of the GenoPred repository download `qctool2` latest versions of `plink1` and `plink2`. If you already have these tools installed, you can also change the paths in `config/config.yaml`. However, if you don't have an up-to-date version, this can lead to bugs :cockroach:, so be careful.

If `workflow/scripts/GenoPred` does not exist after running this step, try to clone the repo manually: `git clone git@github.com:intervene-EU-H2020/GenoPred.git workflow/scripts/GenoPred`. Getting an error? Make sure you have access, and your git is [configured to use ssh](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)!

## :globe_with_meridians: Set up snakemake

While you still have access to the internet, you can [install snakemake using conda/mamba](#step-2-install-snakemake). If you will have access to the internet *and* singularity even while working with sensitive data later on, this will typically be the easiest way to go forward.

The other way is to run snakemake using only the prspipe container. When working with singularity, you can create a local container image which will be able to run the pipeline locally (at least the polygenic scoring and genotype harmonization steps).

```
# this will create a local image called prspipe.sif
singularity pull prspipe.sif docker://rmonti/prspipe:0.0.1
```

After that, you can run snakemake interactively by first starting an [interactive shell inside the container](https://sylabs.io/guides/3.7/user-guide/cli/singularity_shell.html), and then running snakemake:

```
# on the host, in the working directory:
singularity shell -c prspipe.sif

# inside the container:
snakemake --help

# should be the container's R:
which R
RScript --version

# to exit the container
exit
```
> :warning: Note: When running with singularity make sure to clear environment variables related to R such as `R_LIBS`, `R_LIBS_SITE` and `R_LIBS_USER`.

Beware, if the version of R above is *not* `R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"`, your environment variables, `.bashrc` or `.bash_profile` might be interfering with R in the container. Try changing the parameter `Rscript:` inside [`config/config.yaml`](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/config.yaml#L21) to force using the container's R (see the comment in that file).

> Note: If you do **not** have singularity you can also install R-packages manually [install R-packages](#step-3-r-packages-and-other-dependencies). You will then also have to install snakemake.

## Container mounts :file_folder:

:warning: When running the container, you have to "mount" any directories that you want access to, unless they are subdirectories of the main working directory. For example, if you want to access genotype data at `/some/random/path/`, you will have to mount this directory (i.e., make it accessible) inside the container. You can do this by specifying mounts on the command-line when running singularity. The command above becomes `singularity shell -c -B /some/random/ prspipe.sif`. This would make `/some/random` and all its sub-directories available inside the container. If you are running snakemake with the `--singularity` flag (not covered in this tutorial), you have to use the corresponding snakemake parameter, e.g., `--singularity-args "-B /some/random/"`.

## :globe_with_meridians: Download and process the 1000 Genomes data

Assuming you have snakemake running, you can now download and pre-process the 1000 Genomes data. This data will serve as an external reference throughout the analysis and is used for ancestry scoring. First we perform a "quiet" (`-q`) "dry-run" (`-n`).

```
export SNAKEMAKE_CORES=1
bash run.sh -n -q all_setup
```
The command above will display a summary of what will be done. To then actually run these steps:

```
bash run.sh all_setup
```
> :warning: Note: depending on your download speed and how many cores you have available this can take several hours. It will also require a lot of disk-space (~90G).

Once you have successfully completed these steps, you can clear up space by running

```
bash run.sh cleanup_after_setup
```

# Running PRS methods
## Test your setup by running pruning & thresholding on synthetic data

The pipeline ships with synthetic data (genotypes, phenotypes and summary statistics) generated by Sophie Wharrie, Vishnu Raj and Zhiyu Yang. By default, the pipeline is configured to work with this synthetic data. Run the following rule, to extract the synthetic data:

```
bash run.sh initialize_synthetic
```

You can now create polygenic scores using plink and summary statistics from a GWAS run on synthetic data:

```
# first, check what will be executed:
bash run_cluster.sh -n -q prs/pt_clump/synth01/ok

# second, run the rule
bash run_cluster.sh prs/pt_clump/synth01/ok
```

`config/studies.tsv` is the sample-sheet that contains all the information on available summary statistics:

| study_id | ancestry | n_cases | n_controls | ftp_address | local_path | binary | name |
| --- | --- | --- | --- | --- | --- | --- | --- |
| The GWAS acccession | The ancestry abbreviation (EUR, EAS, AMR, SAS and AFR) | The number of case samples | The number of control samples | The GWAS catalog ftp address of the ```.h.tsv.gz``` harmonized summary statistics file, given in the form ```/pub/databases/gwas/summary_statistics.../harmonised/....h.tsv.gz```) | alternatively a local path to a "munged" summary statistics file | phenotype is binary (yes/no) | name |

[`config/studies.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/studies.tsv) is configured to work with summary statistics generated from synthetic data. In order to create scores with a specific PRS method "`{method}`" and GWAS stummary statistics `"{study}"`, the user can request output files which follow the pattern: `prs/{method}/{study}/ok`. Available methods are `dbslmm`,`lassosum`,`ldpred2`,`megaprs`,`prscs`,`pt_clump`.

>:warning:Note: Collaborators don't have to run all these methods. **Their outputs will be distributed.** See steps below.

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


# :globe_with_meridians: Everything below is under (re-)construction :building_construction:. Ignore! 

~~### Download Pre-adjusted summary statistics~~

I've generated adjusted summary statistics for 5 phenoypes (BMI, T2D, breast cancer, prostate cancer and HbA1c). Follow the steps below to download them, and set up the pipeline to use them. 

1.  :globe_with_meridians: Download data from figshare by running `bash run.sh download_test_data`. This might take a while...
2.  Verify the data is in the correct location by running `bash run.sh -n validate_setup_ext`. You should see a message that says *"Nothing to be done."*.

All the steps that require internet access are done. 

~~### Set up Genotype and Phenotype data~~

These steps have not yet been automated, but we will work on automating them in the future. Replace "{bbid}" (biobank id) with a name of your choice in the steps below. **The harmonization of genetic data relies on rsIDs. If your genotypes are not annotated with rsIDs, you will not be able to harmonize your data!**. Contact me if this is the case.

1.  Create folders `custom_input/genotypes/{bbid}` and `custom_input/phenotypes/{bbid}`
2.  :rotating_light: The pipeline requires genotypes in plink format. If you are starting from other formats such as BGEN or VCF, you will first have to convert your data to plink format using [plink v1.9](https://www.cog-genomics.org/plink/1.9/input#oxford) or [QCTOOL v2](https://www.well.ox.ac.uk/~gav/qctool_v2/index.html). If your plink-formated data is not split by chromosome, you can split the data by chromosome with the following loop on the command-line, assuming your genetic data is called `master.{bed,bim,fam}`:
```
# replace "master" with the prefix of your data
# this will generate files tmp_chr{1-22}.{bed,bim,fam}
for i in {1..22}; do plink --bfile master --chr $i --make-bed --out tmp_chr${i}; done
```
3.    :rotating_light: Harmonize the per-chromosome genotype data with the HapMap3 (hm3) variants. The script [hm3_harmoniser.R](https://github.com/intervene-EU-H2020/GenoPred/blob/1d5fddc6e6bf41c7ee94041f84ac91c1afd694fb/Scripts/hm3_harmoniser/hm3_harmoniser.R) can be used to harmonize plink-formatted genotype files with the hm3 variants. If you completed the [basic setup steps above](#basic-setup) successfully, the folder `resources/Geno_1KG` should contain harmonized genotype files with prefixes `1KGPhase3.w_hm3.chr...`. These can be used as the input for the `hm3_harmoniser.R` script (see the `--ref` parameter). Set the `--target` parameter to the prefix of your per-chromosome genotype data prefix (for example `tmp_chr` if you used the loop above to generate them).
> Note: If your genotype panel does not cover >90% of the hm3 variants, imputation might be required!
4.  :rotating_light: Place the per-chromosome harmonized plink-formated bim/bed/fam-files generated in step 3 in  `custom_input/genotypes/{bbid}/`, name them `chr1.{bed,bim,fam}`,`chr2.{bed,bim,fam}`,...,`chr22.{bed,bim,fam}`.
5.  :rotating_light: Extract the phenotypes of interest, in this case BMI, HbA1c, prostate cancer, breast cancer and Type 2 Diabetes. This will be highly specific to the biobank you are working with and may include self-reported endpoints or data scraped from patient records. Phenotypes should be placed in separate files with three tab-separated columns (no header): The family ID, the individual ID, and the Phenotype value (see the plink "pheno" format [here](https://www.cog-genomics.org/plink/1.9/input#pheno)), for example

    |  |  |  |
    | --- | --- | --- |
    | fid_1 | iid_1 | 0 |
    | fid_2 | iid_2 | 1 |
    | ... | ... | ... |

    The IDs should have matches in the plink `.fam`-files generated above, but it is **not** necessary to have phenotype values for all genotyped individuals. This   way only a sub-set of individuals can be analysed for each phenotype to save computational time. Also, if you do not want to analyse a specific phenotype, you can delete the row corresponding to that phenotype from [`conf/studies.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/6b4611bdd2417072229762444a994e02dba6c597/config/studies.tsv).

    Place the phenotype files in `custom_input/phenotypes/{bbid}/`. Name them `{phenotype}.txt`, where {phenotype} should match the entries in the "name"-column of [`config/studies.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/6b4611bdd2417072229762444a994e02dba6c597/config/studies.tsv), i.e `./custom_input/phenotypes/{bbid}/HbA1c.txt`, `./custom_input/phenotypes/{bbid}/BMI.txt` and so on.

    The pipeline does not need to be configured to look for files in `custom_input/phenotypes/*/` and `custom_input/genotypes/*/`. `config/studies.tsv` can be left unchanged if you wish to analyse all phenotypes. Otherwise delete the rows corresponding to the phenotypes you do not wish to analyse from [`conf/studies.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/6b4611bdd2417072229762444a994e02dba6c597/config/studies.tsv).

6.  :rotating_light: Assuming you have downloaded pre-adjusted summary statistics as described above (`bash run.sh download_test_data`), you can now perform hyper-parameter tuning (model selection) on your data. First, perform a dryrun to check that the right jobs will be executed:
    
    ```
    bash run.sh --use-singularity --dryrun --keep-going all_get_best_models_ext
    ```

    the output should look something like this

    ```
    [...]
    Job counts:
            count   jobs
            1       all_get_best_models_ext
            1       ancestry_scoring_ext
            22      calculate_maf_ancestry_ext
            5       dbslmm_score_ext_ref1kg
            5       get_best_models_ext
            5       lassosum_score_ext
            5       ldpred2_score_ext_refukbb
            5       ldpred_score_ext_ref1kg
            5       model_eval_ext
            5       model_eval_ext_prep
            5       prscs_score_ext_refukbb
            5       sbayesr_score_ext_refukbb_robust
            5       sblup_score_ext_ref1kg
            5       sparse_thresholding_score_ext_ref1kg
            79
    ```
    if you see significantly more jobs listed, your pipeline is misconfigured. Please contact me in this case (remo.monti@hpi.de).

This will run ancestry scoring, identify individuals with EUR ancestry, and perform predictions and hyper-parameter tuning for those individuals. You can now run these steps with `bash run.sh --use-singularity --keep-going all_get_best_models_ext`.


    
# Documentation

This part of the documentation will lead you through all steps of the pipeline, including setting it up for new phenotypes.

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

     conda create -c bioconda -c conda-forge -n snakemake snakemake
> Note: make sure to update conda, or use mamba to get the latest version of snakemake i.e. `6.2.1`

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Run `download_resources.sh` and `install_software.sh` to set up software.

    bash ./install_software.sh

### Step 3: R packages and other dependencies

The pipeline relies heavily on [R](https://www.r-project.org/). There are two options to handle R-package dependencies: (1, **recommended**) running the pipeline using singularity containers, or (2) running with a local conda and R installation and installing R-packages manually.

If manually installing R-packages, the following commands will install all dependencies:

    R
    install.packages(c('data.table','doMC','optparse','foreach','caret','ggplot2','cowplot','glmnet','MLmetrics','e1071','stringr','verification', 'RcppArmadillo', 'Matrix', 'fdrtool', 'psych', 'bigsnpr', 'bigreadr', 'runonce'), dependencies=TRUE)
    
    install.packages("bin/lassosum/lassosum_0.4.5.tar.gz", repos=NULL, type="source")

#### Containers (Docker/Singularity)

We have provided [Singularity](https://sylabs.io/) and [Docker](https://www.docker.com/) container [definitions](https://github.com/intervene-EU-H2020/prspipe/tree/main/containers). Users can either build the images from scratch (using `docker build` or `singularity build`), or use the available docker image on [dockerhub](https://hub.docker.com/repository/docker/rmonti/prspipe).

Creating a singularity image that is able to run the pipeline is as simple as:

```
# build from the local container definition...
singularity build containers/singularity/prspipe.sif containers/singularity/prspipe_alldeps.def

# ... or build directly from dockerhub
singularity build containers/singularity/prspipe.sif docker://rmonti/prspipe:0.0.1
```

the latter should work even *without* root access.

To run the pipeline with singulartiy, use the `--use-singularity` flag with snakemake. The default image is defined in `workflow/Snakefile`. Currently, the default is `docker://rmonti/prspipe:0.0.1`. When running with singularity, make sure to clear environment variables `R_LIBS`, `R_LIBS_SITE` and `R_LIBS_USER`, if set, as they can interfere with R in the container.

Instead of running with `--use-singularity` the entire pipeline can also be run locally within the `docker://rmonti/prspipe:0.0.1`-container, as it contains both conda and snakemake.

For example:

```
# build container directly from dockerhub
singularity build containers/singularity/prspipe.sif docker://rmonti/prspipe:0.0.1

# run pipeline commands inside the container
# for example, a dry-run of the setup rules:
singularity run containers/singularity/prspipe.sif ./run.sh -n all_setup
```

### Step 4: Configure workflow

Basic parameters and paths are configured within `config/config.yaml`. The most important input files are described below.

### Step 5: Execute workflow

Activate the snakemake conda environment, and use `./run.sh` to launch snakemake with default parameters:

    conda activate snakemake
    # run all the setup rules (steps 1-3 of GenoPred)
    ./run.sh all_setup
    
    # OR, using singularity
    ./run.sh --use-singularity all_setup

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details (cluster execution, dependency handling etc.).

## Incorporating new GWAS summary statistics

The code implements rules for automatically downloading and pre-processing summary statistics from the [GWAS catalog](https://www.ebi.ac.uk/gwas/home).

To configure and run the code for downloading summary statistics:
1. Select the studies from the GWAS catalog that you want to use
2. For each study, enter the details in the config/studies.tsv file:

| study_id | ancestry | n_cases | n_controls | ftp_address | local_path | binary | name |
| --- | --- | --- | --- | --- | --- | --- | --- |
| The study acccession number given in the GWAS catalog | The ancestry abbreviation (we currently support EUR, EAS, AMR, SAS and AFR) | The number of case samples | The number of control samples | The ftp address of the ```.h.tsv.gz``` harmonized summary statistics file, given in the form ```/pub/databases/gwas/summary_statistics.../harmonised/....h.tsv.gz```) | alternatively a local path to a "munged" summary statistics file | phenotype is binary (yes/no) | name |
| e.g. GCST000998 | e.g. EUR | e.g. 22233 | e.g. 64762 | e.g. /pub/databases/gwas/summary_statistics/GCST000001-GCST001000/GCST000998/harmonised/21378990-GCST000998-EFO_0001645.h.tsv.gz | e.g. ./munged_ss.tsv.gz | e.g. yes | e.g. CAD |

If summary statistics are not available in the harmonized format, consider using [this script](https://github.com/intervene-EU-H2020/prspipe/blob/548640564d81196308cf411b60838d1a6f8a01d6/workflow/scripts/R/munge_sumstats.R) to convert them to munged format.

3. Run the ```all_QC``` rule in the ```base_sumstats.smk``` rule file. Snakemake will download the summary statistics and run scripts to format and QC the data 

```
./run.sh all_QC
```

## Disease prevalence estimates
Model evaluation requires an estimate of disease prevalence in the general population for dichotomous traits. These can be found in the [`pop_prevalence.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/main/config/pop_prevalence.tsv). Make sure to add a corresponding row in that file when adding new endpoints. 

