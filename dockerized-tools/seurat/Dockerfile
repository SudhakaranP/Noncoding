# Seurat

# start with ubuntu
FROM ubuntu

# install required ubuntu packages
RUN rm -rf /var/lib/apt/lists/* && \
                    apt-get clean && \
                    echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list && \
                    gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && \
                    gpg -a --export E084DAB9 | apt-key add - && \
                    apt-get update -y && \
                    apt-get install -y --fix-missing \
                                        curl \
                                        libcairo2-dev \
                                        libcurl4-openssl-dev \
                                        libxml2-dev \
                                        libgtk2.0-dev \
                                        libssl-dev \
                                        xvfb \
                                        xauth \
                                        xfonts-base \
                                        libxt-dev \
                                        r-base \
                                        r-base-dev && \
                    rm -rf /tmp/downloaded_packages/ /tmp/*.rds && \
                    rm -rf /var/lib/apt/lists/*

# set CRAN repository
RUN echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site

# copy r script which will install required packages
ADD setup.r /tmp

# run r script
RUN Rscript /tmp/setup.r

# create data directory (this will be working directory)
RUN mkdir /data && \
                    mkdir /data/input && \
                    mkdir /data/output

# set working directory
WORKDIR /data/output
