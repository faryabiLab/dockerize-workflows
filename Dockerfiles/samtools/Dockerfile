#############################################################
# Dockerfile to build a tool container for samtools
#############################################################

# Set the base image to ubuntu
FROM ubuntu:22.04
LABEL maintainer="jamiesheng0306@gmail.com"
WORKDIR /tmp

# Setup packages
RUN apt-get update && apt-get -y install wget make bzip2 autoconf automake gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 -O samtools.tar.bz2
RUN tar -xjvf samtools.tar.bz2
RUN cd samtools-1.18 && make && make prefix=/usr/local/bin install
RUN ln -s /usr/local/bin/bin/samtools /usr/bin/samtools
