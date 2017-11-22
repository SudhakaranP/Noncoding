# TRANSRATE

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
                                        wget && \
                    rm -rf /tmp/downloaded_packages/ /tmp/*.rds && \
    rm -rf /var/lib/apt/lists/*
    
# make /opt the working directory
WORKDIR /opt

# download and extract transrate binaries
RUN wget https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz && \
	tar zxf transrate-1.0.3-linux-x86_64.tar.gz && \
	rm transrate-1.0.3-linux-x86_64.tar.gz

# add transrate to path
ENV PATH "$PATH:/opt/transrate-1.0.3-linux-x86_64"

# create output directory and make it the working directory
RUN mkdir /data
WORKDIR /data
