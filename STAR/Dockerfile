#############################################################
# Dockerfile to build a sample tool container for STAR
#############################################################

# Set the base image to ubuntu
FROM ubuntu:22.04
LABEL maintainer="jamiesheng0306@gmail.com"
WORKDIR /tmp
ENV star_version 2.7.10b
# Setup packages
RUN apt-get update && apt-get -y install wget build-essential zlib1g-dev
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
RUN tar -xzf ${star_version}.tar.gz
RUN cp STAR-${star_version}/bin/Linux_x86_64/* /usr/local/bin
