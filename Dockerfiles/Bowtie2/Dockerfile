#############################################################
# Dockerfile to build a tool container for Bowtie2
#############################################################

# Set the base image to ubuntu
FROM faryabilab/samtools:0.1.0
LABEL maintainer="jamiesheng0306@gmail.com"
WORKDIR /tmp
# Setup packages
RUN apt-get update && apt-get -y install bowtie2


