FROM ubuntu:16.04

RUN apt-get update
RUN apt-get install -y python-minimal python-pip
RUN apt-get install -y git

# Install Opossum

ENV DST /opossum
RUN git clone https://github.com/BSGOxford/Opossum $DST
RUN mv $DST/* .
RUN rm -rf $DST

COPY ./assets/requirements.txt /requirements.txt
RUN pip install -r requirements.txt
