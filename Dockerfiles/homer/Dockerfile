#############################################################
# Dockerfile to build a tool container for HOMER
# (Hypergeometric Optimization of Motif EnRichment)
#############################################################

# Set the base image to Ubuntu
FROM ubuntu:22.04
LABEL maintainer="nburget12@gmail.com"
WORKDIR /tmp

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

#############################################################
# Install system and development dependencies
#############################################################
RUN apt-get update && apt-get -y install \
    wget curl make gcc g++ perl zip unzip gzip bzip2 \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev \
    libssl-dev libncurses5-dev libgd-dev libxml-simple-perl \
    libdbi-perl libdbd-mysql-perl r-base samtools && \
    rm -rf /var/lib/apt/lists/*

#############################################################
# Download and install HOMER from official source
#############################################################
#############################################################
# Download and install HOMER from official source
#############################################################
WORKDIR /opt/homer

RUN wget http://homer.ucsd.edu/homer/configureHomer.pl
RUN printf "\n" | perl configureHomer.pl -install
#############################################################
# Add HOMER to PATH and environment
#############################################################
ENV PATH="/opt/homer/bin:${PATH}"
ENV HOMER_DIR="/opt/homer"

#############################################################
# Optional: install common genome packages (comment out if not needed)
# RUN perl configureHomer.pl -install hg38
# RUN perl configureHomer.pl -install mm10
#############################################################
