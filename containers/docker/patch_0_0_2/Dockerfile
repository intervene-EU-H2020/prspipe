FROM rmonti/prspipe:0.0.1

SHELL ["/bin/bash","-c"]

# Download plink2 non-avx version
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20220121.zip \
&&  unzip plink2_linux_x86_64_20220121.zip && mv plink2 /usr/local/bin/plink2_x86_64 && rm plink2_linux_x86_64_20220121.zip

# Download plink 2 avx version
RUN /usr/bin/wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20220121.zip \
&&  unzip plink2_linux_avx2_20220121.zip && mv plink2 /usr/local/bin/ && rm plink2_linux_avx2_20220121.zip

# Download plink1 
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip \
&&  unzip plink_linux_x86_64_20201019.zip && mv plink /usr/local/bin/plink

# Install qctool 2
RUN  wget https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.6-Ubuntu16.04-x86_64.tgz \
&&   tar -xvf qctool_v2.0.6-Ubuntu16.04-x86_64.tgz && rm qctool_v2.0.6-Ubuntu16.04-x86_64.tgz \
&&   mv qctool_v2.0.6-Ubuntu16.04-x86_64/qctool /usr/local/bin/ && rm -r qctool_v2.0.6-Ubuntu16.04-x86_64


