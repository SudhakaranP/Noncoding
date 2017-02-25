# Cell Ranger

# Cell Ranger requires Illumina bcl2fastq2 2.17 to be installed on the system.
# bcl2fastq2 was developed and tested on centos6 and RHEL5; we'll use centos6:
FROM centos:6

# install unzip (required for unpackaing bcl2fastq2)
RUN yum install -y unzip

# download and install Illumina bcl2fastq2 version 2.17
#WORKDIR /tmp
#RUN curl -O ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/bcl2fastq/bcl2fastq2-v2.17.1.14-Linux-x86_64.zip && \
    #unzip bcl2fastq2-v2.17.1.14-Linux-x86_64.zip && \
    #rm -f bcl2fastq2-v2.17.1.14-Linux-x86_64.zip && \
    #yum install -y bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm && \
    #rm -f bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm


# download and install Cell Ranger
#WORKDIR /opt
#RUN curl -ko cellranger-1.3.1.tar.gz "https://s3-us-west-2.amazonaws.com/10x.downloads/cellranger-1.3.1.tar.gz?AWSAccessKeyId=AKIAJAZONYDS6QUPQVBA&Expires=1488027285&Signature=ECmZCuIwtpns0MeMvndnZFEhMCg%3D" && \
#    tar -xzf cellranger-1.3.1.tar.gz && \
#    rm -f cellranger-1.3.1.tar.gz

ADD cellranger-1.3.1 /opt/cellranger-1.3.1
ADD anaconda-cr-cs /opt/cellranger-1.3.1/anaconda-cr-cs
ADD cellranger-cs /opt/cellranger-1.3.1/cellranger-cs
ADD cellranger-tiny-ref /opt/cellranger-1.3.1/cellranger-tiny-ref

ENV PATH=/opt/cellranger-1.3.1:$PATH

# create a directory for data and set this as working directory
RUN mkdir /data
WORKDIR /data

