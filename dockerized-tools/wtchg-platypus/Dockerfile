
# base image from python 2.7
FROM python:2.7

ENV PYTHONUNBUFFERED 1

RUN mkdir /data

RUN apt-get update
RUN apt-get -y install wget build-essential zlib1g-dev libncurses5-dev unzip python-dev
RUN pip install cython

ENV C_INCLUDE_PATH=/usr/local/include
ENV LIBRARY_PATH=/usr/local/lib
ENV LD_LIBRARY_PATH=/usr/local/lib

RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2
RUN bunzip2 htslib-1.3.2.tar.bz2
RUN tar -xf htslib-1.3.2.tar
WORKDIR /htslib-1.3.2
RUN ./configure
RUN make
RUN make install
WORKDIR /
RUN rm -r /htslib-1.3.2*

RUN wget http://www.well.ox.ac.uk/bioinformatics/Software/Platypus-latest.tgz
RUN tar -zxf Platypus-latest.tgz
RUN rm Platypus-latest.tgz
RUN mv Platypus* Platypus
WORKDIR /Platypus
RUN ./buildPlatypus.sh

WORKDIR /data

#Clean up APT when done.
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean
RUN apt-get autoremove -y 
RUN rm -rf /var/lib/{apt,dpkg,cache,log}/

