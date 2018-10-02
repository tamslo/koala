# Inspired by https://github.com/BioContainers/containers/blob/master/gatk/3.4/Dockerfile

FROM biocontainers/biocontainers:latest

# INSTALL DEPENDENCIES

USER root

RUN apt-get update -qq && \
    apt-get install -y -q --fix-missing \
      libncurses5-dev  \
      libncursesw5-dev \
      build-essential \
      bzip2 \
      zlib1g-dev \
      maven && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# RUN add-apt-repository ppa:openjdk-r/ppa
# RUN apt-get update
# RUN apt-get install openjdk-7-jdk

# INSTALL SAMTOOLS

ENV ZIP=samtools-1.2.tar.bz2
ENV URL=https://github.com/BioDocker/software-archive/releases/download/samtools
ENV FOLDER=samtools-1.2
ENV DST=/tmp

RUN wget $URL/$ZIP -O $DST/$ZIP && \
    tar xvf $DST/$ZIP -C $DST && \
    rm $DST/$ZIP && \
    cd $DST/$FOLDER && \
    make && \
    make install && \
    cd / && \
    rm -rf $DST/$FOLDER

# Add $HOME/bin to path

ENV PATH=$PATH:/home/biodocker/bin
ENV HOME=/home/biodocker

# INSTALL BWA

ENV ZIP=bwakit-0.7.12_x64-linux.tar.bz2
ENV URL=https://github.com/BioDocker/software-archive/releases/download/bwa
ENV FOLDER=bwa.kit
ENV DST=/home/biodocker/bin

RUN wget $URL/$ZIP -O $DST/$ZIP && \
    tar xvf $DST/$ZIP --strip-components=1 -C $DST && \
    rm $DST/$ZIP

# INSTALL PICARD

ENV ZIP=picard-tools-1.139.zip
ENV URL=https://github.com/broadinstitute/picard/releases/download/1.139/
ENV FOLDER=picard-tools-1.139
ENV DST=/home/biodocker/bin/

RUN wget $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    cd $DST/$FOLDER && \
    mv * .. && \
    cd / && \
    bash -c 'echo -e "#!/bin/bash\njava -jar '$DST'/picard.jar \$@" > '$DST'/picard' && \
    chmod +x $DST/picard && \
    rm -rf $DST/$FOLDER

# INSTALL GATK

ENV ZIP=gatk-4.0.9.0.zip
ENV URL=https://github.com/broadinstitute/gatk/releases/download/4.0.9.0/
ENV FOLDER=gatk-4.0.9.0
ENV DST=/home/biodocker/bin

RUN wget $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    cd $DST/$FOLDER && \
    mv * .. && \
    cd / && \
    chmod +x $DST/gatk && \
    rm -rf $DST/$FOLDER
