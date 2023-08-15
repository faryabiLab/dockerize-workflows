#############################################################
# Dockerfile to build a tool container for Subread
#############################################################

# Set the base image to ubuntu
FROM ubuntu:22.04
LABEL maintainer="jamiesheng0306@gmail.com"
WORKDIR /tmp
# Setup packages
RUN apt-get update && apt-get -y install wget
RUN wget https://github.com/ShiLab-Bioinformatics/subread/releases/download/2.0.2/subread-2.0.2-Linux-x86_64.tar.gz
RUN tar xvzf subread-2.0.2-Linux-x86_64.tar.gz && rm subread-2.0.2-Linux-x86_64.tar.gz \
	&& ln -sf /tmp/subread-2.0.2-Linux-x86_64/bin/* /usr/local/bin/
