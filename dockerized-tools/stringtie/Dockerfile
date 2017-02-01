# Stringtie

FROM ubuntu

# create output directory
RUN mkdir /data

# download latest package lists
RUN apt-get update
# install tools which may be needed
RUN apt-get -y install wget build-essential zlib1g-dev libncurses5-dev unzip

RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.2d.Linux_x86_64.tar.gz

RUN tar -zxf stringtie-1.3.2d.Linux_x86_64.tar.gz

RUN cp /stringtie-1.3.2d.Linux_x86_64/stringtie /usr/local/bin/

RUN rm -r /stringtie-1.3.2d.Linux_x86_64*

WORKDIR /data


