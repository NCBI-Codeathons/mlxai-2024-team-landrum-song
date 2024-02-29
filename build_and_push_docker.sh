#!/bin/bash

# build the docker image with a tag
docker build -t nrminor/ncbi-ml-ai:v0.0.1 .

# push to docker hug
docker push nrminor/ncbi-ml-ai:v0.0.1

# example docker run command
# docker run \
# --user $(id -u):$(id -g) -it -v $(pwd):/scratch \
# -w /scratch nrminor/ncbi-ml-ai:v0.0.1
