# Basically from https://github.com/Illumina/hap.py/blob/d51d111e494b561b37c66299daf5a6c65a8d2ca9/Dockerfile
FROM ubuntu:14.04

RUN apt-get update
RUN apt-get upgrade -y

RUN apt-get install -y python3-minimal git
RUN apt-get install -y python2.7 python2.7-dev python-pip cython
RUN apt-get install -y python-setuptools python-numpy zlib1g-dev
RUN easy_install -U distribute
RUN apt-get install -y build-essential libatlas-base-dev gfortran pkg-config
RUN apt-get install -y software-properties-common python-software-properties
RUN apt-get install -y cmake libncurses5-dev libbz2-dev liblzma-dev
RUN apt-get install -y autoconf wget

RUN pip install --upgrade cython
RUN pip install --upgrade numpy
RUN pip install --upgrade pandas
RUN pip install pybedtools
RUN pip install pysam
RUN pip install bx-python
RUN pip install scipy

RUN echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
    add-apt-repository -y ppa:webupd8team/java && \
    apt-get update && \
    apt-get install -y oracle-java8-installer && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/cache/oracle-jdk8-installer

WORKDIR /opt
RUN wget http://archive.apache.org/dist/ant/binaries/apache-ant-1.9.7-bin.tar.gz && \
    tar xzf apache-ant-1.9.7-bin.tar.gz && \
    rm apache-ant-1.9.7-bin.tar.gz
ENV PATH $PATH:/opt/apache-ant-1.9.7/bin
WORKDIR /

# Install Haplotype Comparison Tools

ENV DST /opt/hap.py-source
RUN git clone https://github.com/Illumina/hap.py.git $DST

WORKDIR $DST
ENV DST /opt/hap.py
RUN python install.py $DST --with-rtgtools --no-test
WORKDIR /

# Download GIAB high confidence calls

COPY ./assets/download_calls.py ./download_calls.py
RUN python3 download_calls.py

COPY ./assets/evaluate_variants.sh ./evaluate_variants.sh
