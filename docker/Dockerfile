FROM debian:stable

LABEL maintainer "Philipp Kirstahler" \
      version "0.1" \
      description "plasmidProject"

# install needed tools
RUN apt-get update --fix-missing -qq && apt-get install -y -q \
    curl \
    locales \
    libncurses5-dev  \
    libncursesw5-dev \
    build-essential \
    pkg-config \
    zlib1g-dev \
    bzip2 \
    r-base \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


RUN apt-get update && apt-get install -y python

# install SAMtools
RUN curl -fksSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar xj && \
    cd samtools-1.3.1 && \
    make all all-htslib && make install install-htslib

# install bedtools
RUN curl -fksSL https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz | tar xz && \
    cd bedtools2 && \
    make && make install

# install minimap2
RUN git clone https://github.com/lh3/minimap2.git && \
    cd minimap2 && \
    make

ENV PATH "$PATH:/minimap2"


# install minimasm
RUN git clone https://github.com/lh3/miniasm.git && \
    cd miniasm && \
    make

ENV PATH "$PATH:/miniasm"

## Install R packages for tidyverse
RUN R -e 'install.packages( c("tidyverse"), repos="http://cloud.r-project.org/");'

