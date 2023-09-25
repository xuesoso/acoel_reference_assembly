#!/bin/bash

# Build the docker image with the tag "2023_acoel_ref" and the same user name and group permission as the host user
docker build Docker -f Docker/Dockerfile -t 2023_acoel_ref --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g)
