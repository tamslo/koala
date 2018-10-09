# Inspired by https://github.com/AveraSD/ngs-docker-star
FROM ubuntu:16.04

ENV STAR_VERSION 2.6.0c

# From https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
RUN apt-get update
RUN apt-get install -y g++
RUN apt-get install make

# Other required packages
RUN apt-get install -y wget
RUN apt-get install zlib1g-dev

RUN wget -P /opt/star https://github.com/alexdobin/STAR/archive/$STAR_VERSION.tar.gz
RUN tar -xzf /opt/star/$STAR_VERSION.tar.gz -C /opt/star
RUN make STAR -C /opt/star/STAR-$STAR_VERSION/source

#Add star to PATH
ENV PATH /opt/star/STAR-$STAR_VERSION/bin/Linux_x86_64:$PATH
