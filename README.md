# Snakemake workflow: prspipe

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/prspipe.svg?branch=master)](https://travis-ci.org/snakemake-workflows/prspipe)

Snakemake pipeline to run Polygenic Risk Score (PRS) prediction. Implements and extends the [GenoPred](https://github.com/opain/GenoPred) pipeline, i.e. a reference standardized framework for the prediction of PRS using different state of the art methods using summary statistics.

## Authors

* Remo Monti & Sophie Wharrie
* [GenoPred](https://github.com/opain/GenoPred): Oliver Pain 

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

## Quick Start for Collaborators

The full pipeline can roughly be devided into two stages: (1) Download and adjustment of summary statistics using pre-computed LD reference panels where available and data from the 1000 Genomes project (publica data), and (2) prediction and evaluation of PRS using the adjusted summary statistics, which includes hyper-parameter tuning using cross validation. The second step uses private data e.g. from the UK Biobank or your biobank of interest. **We are currently evaluating running the second step in other biobanks**.

Rules that re-implement the analysis of the UK Biobank data as shown in the GenoPred paper can be found in [`workflow/rules/ukbb.smk`](https://github.com/intervene-EU-H2020/prspipe/blob/main/workflow/rules/ukbb.smk). Follow the steps below in order to work with data from other biobanks using rules defined [`workflow/rules/external_biobanks.smk`](https://github.com/intervene-EU-H2020/prspipe/blob/main/workflow/rules/external_biobanks.smk).

Steps that need internet access are marked with :globe_with_meridians: and steps that require access to sensitive data are marked with :rotating_light:.

### Basic Setup

1.  :globe_with_meridians: Install conda and [snakemake](#step-2-install-snakemake), if not available already.
2.  :globe_with_meridians: Clone the repository, and run `bash ./install_software.sh`. This will download necessary software dependencies. Make sure you have access to https://github.com/intervene-EU-H2020/GenoPred and your git is configured with ssh, otherwise this will not work!
    1.  If `workflow/scripts/GenoPred` does not exist after running this step, try to clone the repo manually: `git clone git@github.com:intervene-EU-H2020/GenoPred.git workflow/scripts/GenoPred`. Getting an error? Make sure you have access, and your git is [configured to use ssh](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)!
4.  :globe_with_meridians: If singularity is **not** available, [install R-packages](#step-3-r-packages-and-other-dependencies).
5.  :globe_with_meridians: Download resources by running `bash run.sh --use-singularity get_plink_files_chr_all download_hapmap3_snplist`.
6.  Process the 1000 Genomes data by running `bash run.sh --use-singularity all_setup`
    > Careful, this will need quite some disk space (almost 90G)
7.  If step 6 completed successfully, you can clear up space by running `bash run.sh cleanup_after_setup`, which will remove many intermediate files.

When running with singularity, make sure to clear environment variables related to R such as `R_LIBS`, `R_LIBS_SITE` and `R_LIBS_USER`. Step 5 should be run on a compute-node.

### Download Pre-adjusted summary statistics

I've generated adjusted summary statistics for 5 phenoypes (BMI, T2D, breast cancer, prostate cancer and HbA1c). Follow the steps below to download them, and set up the pipeline to use them. 

1.  :globe_with_meridians: Download data from figshare by running `bash run.sh download_test_data`. This might take a while...
2.  Verify the data is in the correct location by running `bash run.sh -n validate_setup_ext`. You should see a message that says *"Nothing to be done."*.

All the steps that require internet access are done. 

### Set up Genotype and Phenotype data

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

## Downloading GWAS summary statistics

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
