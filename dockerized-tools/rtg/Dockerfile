# Create image for RTG Tools
# base image from openjdk 8
FROM openjdk:8

# create output directory
RUN mkdir /data
# create directory to hold sdf files
RUN mkdir /sdf

# download latest package lists
RUN apt-get update
# install tools which may be needed
RUN apt-get -y install wget build-essential zlib1g-dev libncurses5-dev unzip ant

# get RTG tools source code and build
RUN git clone https://github.com/RealTimeGenomics/rtg-tools.git
WORKDIR /rtg-tools
RUN ant runalltests
RUN ant zip-nojre
# install RTG tools
WORKDIR /opt
RUN unzip /rtg-tools/dist/rtg-tools-3.7.1-eb13bbb-nojre.zip
RUN mv rtg-tools-* rtg-tools
RUN rm -r /rtg-tools

# copy configuration file to installation directory - this prevents RTG from prompting user with question on first launch
ADD rtg.cfg /opt/rtg-tools

# set working directory
WORKDIR /data
