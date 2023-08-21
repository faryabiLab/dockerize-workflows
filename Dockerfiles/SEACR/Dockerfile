#############################################################
# Dockerfile to build a tool container for SEACR
#############################################################

# Set the base image to ubuntu
FROM ubuntu:22.04
LABEL maintainer="jamiesheng0306@gmail.com"

WORKDIR /tmp
# set timezone
RUN ln -snf /usr/share/zoneinfo/$CONTAINER_TIMEZONE /etc/localtime && echo $CONTAINER_TIMEZONE > /etc/timezone
# Setup packages
RUN apt-get update && apt-get -y install wget bedtools r-base r-base-dev
RUN wget https://github.com/FredHutch/SEACR/archive/refs/tags/v1.3.tar.gz -O seacr.tar.gz && \
	tar zxvf seacr.tar.gz && rm seacr.tar.gz

