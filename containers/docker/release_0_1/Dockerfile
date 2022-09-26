
FROM bioconductor/bioconductor_docker:RELEASE_3_14

CMD /bin/bash "$@" 

SHELL ["/bin/bash","-c"]

# Installing miniconda
RUN  touch /.condarc \
&&   wget -O Miniconda_installer.sh -c https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh \
&&   /bin/bash Miniconda_installer.sh -bfp /usr/local \
&&   conda config --file /.condarc --add channels defaults \
&&   conda config --file /.condarc --add channels bioconda \
&&   conda config --file /.condarc --add channels conda-forge

# Installing conda dependancies
RUN  conda install -c conda-forge -c bioconda  bcftools mamba jupyterlab \
&&   mamba install -c conda-forge -c bioconda snakemake=6.2.1 \
&&   conda clean --all

# Download plink2 
RUN wget -O plink2.zip https://figshare.com/ndownloader/files/33895304?private_link=3641ef6df51eddbeea60 \
&&  unzip plink2.zip && mv plink2 /usr/local/bin/plink2 && mv plink2_x86_64 /usr/local/bin/plink2_x86_64 && rm -f plink2.zip

# Download plink1 
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip \
&&  unzip plink_linux_x86_64_20201019.zip && mv plink /usr/local/bin/plink && rm -f plink_linux_x86_64_20201019.zip

# Download qctool 2
RUN  wget https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.6-Ubuntu16.04-x86_64.tgz \
&&   tar -xvf qctool_v2.0.6-Ubuntu16.04-x86_64.tgz && rm qctool_v2.0.6-Ubuntu16.04-x86_64.tgz \
&&   mv qctool_v2.0.6-Ubuntu16.04-x86_64/qctool /usr/local/bin/ && rm -r qctool_v2.0.6-Ubuntu16.04-x86_64

# Installing R packages (1)
# these are the vasic dependencies we always want to be there

RUN R -e "install.packages(c('data.table','doMC','optparse','foreach','ggplot2','cowplot','MLmetrics','stringr','verification', 'RcppArmadillo', 'Matrix', 'fdrtool', 'psych', 'bigsnpr', 'bigreadr', 'runonce', 'GGally', 'NbClust', 'IRkernel'))"

# enable jupyter R kernel
RUN R -e "library(IRkernel); IRkernel::installspec(name = 'r', displayname = 'R')"

RUN wget https://github.com/tshmak/lassosum/releases/download/v0.4.5/lassosum_0.4.5.tar.gz \
&&   R -e "install.packages('lassosum_0.4.5.tar.gz', repos=NULL, type='source')" \
&&   rm lassosum_0.4.5.tar.gz

# Installing R packages (2)
# caret needs to be installed with dependencies

RUN  R -e "install.packages('caret', dependencies=TRUE)"

USER rstudio

# make the R-kernel available in jupyterlab for user rstudio
RUN  R -e "library(IRkernel); IRkernel::installspec(name = 'r', displayname = 'R')"

USER root

# Installing R packages (3)
# install machine learning packages, easily added as additional layers later...

RUN R -e "install.packages(c('glmnet', 'e1071', 'RSQLite'))"


