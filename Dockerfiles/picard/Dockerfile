#############################################################
# Dockerfile to build a tool container for picard
#############################################################

# Set the base image to samtools
FROM faryabilab/samtools:0.1.0
LABEL maintainer="jamiesheng0306@gmail.com"
WORKDIR /tmp
ENV version 3.1.0
ENV PICARD /usr/local/bin/picard/picard.jar
WORKDIR /usr/local/bin/picard/

RUN apt-get update -y && apt-get install -y \
    openjdk-17-jre \
    curl \
    less \
    vim  \
    wget

RUN mkdir -p /usr/local/bin/picard \
    && curl -SL https://github.com/broadinstitute/picard/releases/download/${version}/picard.jar \
    > /usr/local/bin/picard/picard.jar

WORKDIR /tmp
