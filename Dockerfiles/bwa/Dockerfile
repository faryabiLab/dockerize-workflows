#############################################################
# Dockerfile to build a tool container for bwa
#############################################################

# Set the base image to ubuntu
FROM ubuntu:22.04
LABEL maintainer="jamiesheng0306@gmail.com"
WORKDIR /tmp
# Setup packages
RUN apt-get update && apt-get -y install make zlib1g-dev git build-essential
RUN git clone https://github.com/lh3/bwa.git && cd bwa && make && cp bwa /usr/local/bin

