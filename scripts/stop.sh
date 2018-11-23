#!/bin/bash

killall python
docker stop $(docker ps -aq)
docker rm $(docker ps -aq)
