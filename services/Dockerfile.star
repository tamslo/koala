FROM ubuntu:16.04
WORKDIR /home
RUN echo "Star up and running"
RUN ls
CMD tail -F /dev/null # leave container running
