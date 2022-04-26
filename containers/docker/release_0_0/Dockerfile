
FROM bioconductor/bioconductor_docker:RELEASE_3_13

RUN unzip test

# Installing miniconda
RUN    touch /.condarc \
&&     wget -O Miniconda_installer.sh -c https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh \
&&     /bin/bash Miniconda_installer.sh -bfp /usr/local \
&&     conda config --file /.condarc --add channels defaults \
&&     conda config --file /.condarc --add channels bioconda \
&&     conda config --file /.condarc --add channels conda-forge

# Installing conda dependancies
RUN  conda install -c bioconda bcftools \
&&   conda install -c conda-forge mamba \
&&   mamba install -c conda-forge -c bioconda snakemake=6.2.1 \
&&   conda clean --all

# Installing R packages
RUN  wget https://github.com/tshmak/lassosum/releases/download/v0.4.5/lassosum_0.4.5.tar.gz \ 
&&   R -e "install.packages(c('data.table','doMC','optparse','foreach','caret','ggplot2','cowplot','glmnet','MLmetrics','e1071','stringr','verification', 'RcppArmadillo', 'Matrix', 'fdrtool', 'psych', 'bigsnpr', 'bigreadr', 'runonce'), dependencies=TRUE)" \
&&   R -e "install.packages('lassosum_0.4.5.tar.gz', repos=NULL, type='source')" \
&&   R -e "BiocManager::install('MungeSumstats')" \
&&   rm lassosum_0.4.5.tar.gz

CMD /bin/bash "$@"

SHELL ["/bin/bash"]

