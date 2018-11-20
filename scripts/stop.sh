#!/bin/bash

killall python
docker stop $(docker ps -aq)
docker remove $(docker ps -aq
