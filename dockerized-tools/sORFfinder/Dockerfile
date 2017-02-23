# sORFfinder

# start with ubuntu
FROM ubuntu

# install required ubuntu packages - n.b. for original sORFfinder need: blast2
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get clean
RUN apt-get update -y && apt-get install -y --fix-missing\
                                        wget\
                                        clustalw\
                                        ncbi-blast+
                    

# cleanup
RUN rm -rf /tmp/downloaded_packages/ /tmp/*.rds
RUN rm -rf /var/lib/apt/lists/*

RUN mkdir /sORFfinder2
ADD src/* /sORFfinder2/
RUN chmod 755 /sORFfinder2/*.pl

RUN mkdir /data
WORKDIR /data

