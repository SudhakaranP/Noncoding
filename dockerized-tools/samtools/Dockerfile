#samtools

FROM ubuntu

# https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2

# install required ubuntu packages
RUN rm -rf /var/lib/apt/lists/* && \
     apt-get clean && \
     apt-get update -y && \
     apt-get install -y --fix-missing \
          build-essential \
          zlib1g-dev \
          libncurses5-dev \
          libbz2-dev \
          liblzma-dev \
          wget && \
     rm -rf /tmp/downloaded_packages/ /tmp/*.rds && \
     rm -rf /var/lib/apt/lists/*

# install samtools
WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2 && \
     tar xf samtools-1.4.1.tar.bz2 && \
     rm samtools-1.4.1.tar.bz2
WORKDIR /tmp/samtools-1.4.1
RUN make && \
     make prefix=/usr/local install

# create output directory and make it the working directory
RUN mkdir /data
WORKDIR /data
