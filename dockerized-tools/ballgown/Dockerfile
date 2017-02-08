# Example R tool

# get specific version of R (version 3.3.2 was latest when this file was created)
FROM r-base:3.3.2
# alternatively request latest
#FROM r-base:latest

# change CMD from default of "R" to bash
CMD ["bash"]

# install required ubuntu packages
RUN apt-get update -y
RUN apt-get install curl libcairo2-dev libcurl4-openssl-dev libxml2-dev libgtk2.0-dev xvfb xauth xfonts-base libxt-dev -y

# copy r script which will install required packages
ADD setup.r /tmp

# make r script executable
RUN chmod 755 /tmp/setup.r

# run r script
RUN /tmp/setup.r

# create data directory (this will be working directory)
RUN mkdir /data
RUN mkdir /data/read-coverage-tables
RUN mkdir /data/output

# set working directory
WORKDIR /data



