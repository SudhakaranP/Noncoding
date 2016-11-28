
# base image from python 2.7
FROM python:2.7

ENV PYTHONUNBUFFERED 1

RUN mkdir /data
RUN mkdir /idx
RUN apt-get update
RUN apt-get -y install wget build-essential zlib1g-dev libncurses5-dev unzip

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.5-Linux_x86_64.zip
RUN unzip hisat2-2.0.5-Linux_x86_64.zip
RUN cp hisat2-2.0.5/hisat2* hisat2-2.0.5/*.py /usr/local/bin/
RUN rm hisat2-2.0.5-Linux_x86_64.zip

RUN wget http://downloads.sourceforge.net/project/samtools/samtools/1.3.1/samtools-1.3.1.tar.bz2
RUN tar jxvf samtools-1.3.1.tar.bz2
RUN rm samtools-1.3.1.tar.bz2
WORKDIR /samtools-1.3.1
RUN make
RUN make install
WORKDIR /data

ADD align.py /usr/local/bin/

