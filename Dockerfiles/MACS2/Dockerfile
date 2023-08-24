#############################################################
# Dockerfile to build a tool container for MACS2
#############################################################

# Set the base image to ubuntu
FROM ubuntu:22.04
LABEL maintainer="jamiesheng0306@gmail.com"
WORKDIR /tmp
# Setup packages
RUN apt-get update && apt-get -y install python3-pip
RUN pip install macs2

