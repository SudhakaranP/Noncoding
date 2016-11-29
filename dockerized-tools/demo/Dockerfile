
# base image from python 2.7
FROM python:2.7

ENV PYTHONUNBUFFERED 1

# create a data directory
RUN mkdir /data

# make this data directory the working directory
WORKDIR /data

# add the python script to the bin directory for tools not installed using a package manager
ADD demo.py /usr/local/bin/

