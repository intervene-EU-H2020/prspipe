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

Rules that re-implement the analysis of the UK Biobank data as shown in the GenoPred paper can be found in [`workflow/rules/ukbb.smk`](https://github.com/intervene-EU-H2020/prspipe/blob/main/workflow/rules/ukbb.smk). Follow the steps below in order to work with data from other biobanks.

### Basic Setup

1.  Install conda and [snakemake](#step-2-install-snakemake), if not available already.
2.  Clone the repository, and run `bash ./install.sh`. This will download necessary software dependencies.
3.  If singularity is **not** available, [install R-packages](#step-3-r-packages-and-other-dependencies).
4.  Download resources by running `bash run.sh --use-singularity get_plink_files_chr_all download_hapmap3_snplist`.
5.  Process the 1000 Genomes data by running `bash run.sh --use-singularity all_setup`

When running with singularity, make sure to clear environment variables related to R such as `R_LIBS`, `R_LIBS_SITE` and `R_LIBS_USER`. Step 5 should be run on a compute-node.

### Download Pre-adjusted summary statistics

I've generated adjusted summary statistics for 5 phenoypes (BMI, T2D, breast cancer, prostate cancer and HbA1c). Follow the steps below to download them, and set up the pipeline to use them. 

1.  Download data from figshare by running `bash run.sh download_test_data`. This might take a while...
2.  Verify the data is in the correct location by running `bash run.sh -n validate_setup_ext`. You should see a message that says *"Nothing to be done."*.

### Set up Genotype and Phenotype data

These steps have not yet been automated, but we will work on automating them in the future. Replace "{bbid}" with a suitable name in the steps below.

1.  Create folders `custom_input/genotypes/{bbid}` and `custom_input/phenotypes/{bbid}`
3.  The pipeline requires genotypes in plink format. Harmonize your genotype data with the HapMap3 (hm3) variants. The GenoPred script [Harmonisation_of_UKBB.R](https://github.com/intervene-EU-H2020/GenoPred/blob/1d5fddc6e6bf41c7ee94041f84ac91c1afd694fb/Scripts/Harmonisation_of_UKBB/Harmonisation_of_UKBB.R) illustrates these steps for the UK Biobank data. The script [hm3_harmoniser.R](https://github.com/intervene-EU-H2020/GenoPred/blob/1d5fddc6e6bf41c7ee94041f84ac91c1afd694fb/Scripts/hm3_harmoniser/hm3_harmoniser.R) can be used to harmonize plink-formatted genotype files with the hm3 variants. If you completed the [basic setup steps above](#basic-setup) successfully, the folder `resources/Geno_1KG` should contain harmonized genotype files with prefixes `1KGPhase3.w_hm3.chr...`. These can be used as the input for the `hm3_harmoniser.R` script (see the `--ref` parameter).
> Note: If your genotype panel does not cover >90% of the hm3 variants, imputation might be required!
4.  Place per-chromosome harmonized plink-formated bim/bed/fam-files in  `custom_input/genotypes/{bbid}/`, name them `chr1.bed`,`chr2.bed`,...,`chr22.bed`.
5.  Place phenotype files in `custom_input/genotypes/{bbid}/`. Name them `{phenotype}.txt`, where {phenotype} should match the entries in the "name"-column of [`studies.tsv`](https://github.com/intervene-EU-H2020/prspipe/blob/091a9184130a05942840fab6bb3dc5ede59beb6e/config/studies_new.tsv), i.e `HbA1c.txt`, `BMI.txt`, `T2D.txt`, `Prostate_cancer.txt` and `Breast_cancer.txt`. These files should have 3 columns: The family ID, the individual ID, and the Phenotype value (see also [here](https://www.cog-genomics.org/plink/1.9/input#pheno)).

Assuming you have downloaded pre-adjusted summary statistics (`bash run.sh download_test_data`), you can now perform hyper-parameter tuning (model selection) on your data. Contact me (remo.monti@hpi.de) before trying to run the step below.

```
bash run.sh --use-singularity all_get_best_models_ext
```

This will run ancestry scoring, identify individuals with EUR ancestry, and perform predictions and hyper-parameter tuning for those individuals. 

# Documentation
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

This should work even *without* root access.

To run the pipeline with singulartiy, use the `--use-singularity` flag with snakemake. The default image is defined in `workflow/Snakefile`. Currently, the default is `docker://rmonti/prspipe:0.0.1`. When running with singularity, make sure to clear environment variables `R_LIBS`, `R_LIBS_SITE` and `R_LIBS_USER`, if set, as they can interfere with R in the container.

Instead of running with `--use-singularity` the entire pipeline can also be run locally within the `docker://rmonti/prspipe:0.0.1`-container, as it contains both conda and snakemake.

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
