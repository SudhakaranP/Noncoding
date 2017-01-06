
# base image from python 2.7
FROM python:2.7

ENV PYTHONUNBUFFERED 1

RUN apt-get update
RUN apt-get -y install wget unzip

RUN pip install numpy
RUN pip install PyVCF
RUN pip install PyFASTA

# /opt is probably a reasonable location for the smash benchmarking tools
WORKDIR /opt
RUN wget https://github.com/amplab/smash/archive/master.zip
RUN unzip master.zip
RUN rm master.zip

RUN mkdir /data
WORKDIR /data

