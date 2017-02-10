# Ballgown

# start with ubuntu
FROM ubuntu

# install required ubuntu packages
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get clean
RUN apt-get update -y
RUN apt-get install -y --fix-missing\
                    curl \
                    libcairo2-dev \
                    libcurl4-openssl-dev \
                    libxml2-dev \
                    libgtk2.0-dev \
                    xvfb \
                    xauth \
                    xfonts-base \
                    libxt-dev \
                    r-base \
                    r-base-dev

# cleanup
RUN rm -rf /tmp/downloaded_packages/ /tmp/*.rds
RUN rm -rf /var/lib/apt/lists/*

# set CRAN repository
RUN echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site

# copy r script which will install required packages
ADD setup.r /tmp

# run r script
RUN Rscript /tmp/setup.r

# create data directory (this will be working directory)
RUN mkdir /data
RUN mkdir /data/read-coverage-tables
RUN mkdir /data/results

# set working directory
WORKDIR /data/results
