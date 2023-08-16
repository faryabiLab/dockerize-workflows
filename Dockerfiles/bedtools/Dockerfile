#############################################################
# Dockerfile to build a tool container for bedtools
#############################################################

# Set the base image to samtools
FROM faryabilab/samtools:0.1.0
LABEL maintainer="jamiesheng0306@gmail.com"
WORKDIR /tmp
ENV version 2.31.0

RUN apt-get update -y && apt-get install -y \
    wget \
    bc \
    make \
    build-essential \
    zlib1g-dev  \
    liblzma-dev  \
    libbz2-dev  \
    python3 \
    curl

RUN curl -SLo /usr/local/bin/bedGraphToBigWig http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bedGraphToBigWig && \
chmod +x /usr/local/bin/bedGraphToBigWig

RUN wget https://github.com/arq5x/bedtools2/releases/download/v${version}/bedtools-${version}.tar.gz -O bedtools.tar.gz && \
	tar zxvf bedtools.tar.gz && rm bedtools.tar.gz && \
	cd bedtools2 && \
	make && \
  cp ./bin/* /usr/local/bin/

