#vcflib

FROM ubuntu

# install required ubuntu packages
RUN rm -rf /var/lib/apt/lists/* && \
                    apt-get clean && \
                    apt-get update -y && \
                    apt-get install -y --fix-missing \
                                        build-essential \
                                        zlib1g-dev \
                                        libncurses5-dev \
                                        unzip \
                                        git && \
                    rm -rf /tmp/downloaded_packages/ /tmp/*.rds && \
    rm -rf /var/lib/apt/lists/*
    
WORKDIR /opt
RUN git clone --recursive https://github.com/vcflib/vcflib.git
WORKDIR /opt/vcflib
RUN make openmp
ENV PATH "$PATH:/opt/vcflib/bin"

# create output directory and make it the working directory
RUN mkdir /data
WORKDIR /data
