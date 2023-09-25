#!/bin/bash

# Build the docker image first with "bash build_docker.sh"
# This script will run the docker container, mount the current directory, and pass port 8887 to the container
# The host can access the notebook frontend by going to http://localhost:8887
docker run -it --rm -v $(pwd):/tmp/:Z 2023_acoel_ref bash
