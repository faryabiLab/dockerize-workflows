#############################################################
# Dockerfile to build a sample tool container for Trim Galore
#############################################################

# Set the base image to ubuntu
FROM ubuntu:22.04
LABEL maintainer="jamiesheng0306@gmail.com"

# Setup packages
RUN apt-get update && apt-get -y install curl cutadapt fastqc
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o /tmp/trim_galore.tar.gz
RUN tar xvzf /tmp/trim_galore.tar.gz
RUN mv TrimGalore-0.6.10/trim_galore /usr/local/bin/
