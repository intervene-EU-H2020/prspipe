# Snakemake workflow: prspipe

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/prspipe.svg?branch=master)](https://travis-ci.org/snakemake-workflows/prspipe)

Snakemake pipeline to run Polygenic Risk Score (PRS) prediction. Implements and extends the [GenoPred](https://github.com/opain/GenoPred) pipeline, i.e. a reference standardized framework for the prediction of PRS using different state of the art methods.

## Authors

* Remo Monti & Sophie Wharrie
* [GenoPred](https://github.com/opain/GenoPred): Oliver Pain 

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify your sample setup.

> Note: samples.tsv is currently disabled.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

     conda create -c bioconda -c conda-forge -n snakemake snakemake
> Note: make sure to update conda, or use mamba to get the latest version of snakemake i.e. `6.2.1`

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Run `download_resources.sh` and `install_software.sh` to set up software.

    bash ./download_resources.sh
    bash ./install_software.sh


### Step 4: R packages and other dependencies

The pipeline relies heavily on [R](https://www.r-project.org/). There are two options to handle R-package dependencies: (1) running with a local R installation or and installing packages manually or (2) running the pipeline using singularity containers.

If running with a local R installation, the following commands will install all dependencies:

    R
    install.packages(c('data.table','doMC','optparse','foreach','caret','ggplot2','cowplot','glmnet','MLmetrics','e1071','stringr','verification', 'RcppArmadillo', 'Matrix', 'fdrtool', 'psych', 'bigsnpr', 'bigreadr', 'runonce'), dependencies=TRUE)
    
    install.packages("bin/lassosum/lassosum_0.4.5.tar.gz", repos=NULL, type="source")

#### Containers (Docker/Singularity)

Alternatively, we have provided [Singularity](https://sylabs.io/) and [Docker](https://www.docker.com/) container [definitions](https://github.com/intervene-EU-H2020/prspipe/tree/main/containers). Users can either build the images from scratch (using `docker build` or `singularity build`), or use the available docker image on [dockerhub](https://hub.docker.com/repository/docker/rmonti/prspipe). 

Creating a singularity image that is able to run the pipeline is as simple as:

```
singularity build containers/singularity/prspipe.sif containers/singularity/prspipe_alldeps_fromdocker.def
```

> Remo: Docker and Singularity are not fully supported yet (August 4, 2021)


### Step 5: Execute workflow

Activate the conda environment, and use `./run.sh` to launch snakemake with default parameters:

    conda activate snakemake
    # run all the setup rules (steps 1-3 of GenoPred)
    ./run.sh all_setup

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details (cluster execution, dependency handling etc.).

### Step 6: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.
An example (using some trivial test data) can be seen [here](https://cdn.rawgit.com/snakemake-workflows/rna-seq-kallisto-sleuth/master/.test/report.html).

### Step 7: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 8: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/prspipe.git` or `git remote add -f upstream https://github.com/snakemake-workflows/prspipe.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 9: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

## Testing

Test cases are in the subfolder `.test`. They are automatically executed via continuous integration with [Github Actions](https://github.com/features/actions).

## Downloading GWAS summary statistics

The code implements rules for automatically downloading and pre-processing summary statistics from the [GWAS catalog](https://www.ebi.ac.uk/gwas/home).

To configure and run the code for downloading summary statistics:
1. Select the studies from the GWAS catalog that you want to use
2. For each study, enter the details in the config/studies.tsv file:

| study_id | ancestry | n_cases | n_controls | ftp_address | local_path | binary | name |
| --- | --- | --- | --- | --- | --- | --- | --- |
| The study acccession number given in the GWAS catalog | The ancestry abbreviation (we currently support EUR, EAS, AMR, SAS and AFR) | The number of case samples | The number of control samples | The ftp address of the ```.h.tsv.gz``` harmonized summary statistics file, given in the form ```/pub/databases/gwas/summary_statistics.../harmonised/....h.tsv.gz```) | alternatively a local path to a "munged" summary statistics file | phenotype is binary (yes/no) | name |
| e.g. GCST000998 | e.g. EUR | e.g. 22233 | e.g. 64762 | e.g. /pub/databases/gwas/summary_statistics/GCST000001-GCST001000/GCST000998/harmonised/21378990-GCST000998-EFO_0001645.h.tsv.gz | e.g. ./munged_ss.tsv.gz | e.g. yes | e.g. CAD_Schunkert |

If summary statistics are not available in the harmonized format, consider using [this script](https://github.com/intervene-EU-H2020/prspipe/blob/548640564d81196308cf411b60838d1a6f8a01d6/workflow/scripts/R/munge_sumstats.R) to convert them to munged format.

3. Run the ```all_QC``` rule in the ```base_sumstats.smk``` rule file. Snakemake will download the summary statistics and run scripts to format and QC the data 

```
./run.sh all_QC
```
