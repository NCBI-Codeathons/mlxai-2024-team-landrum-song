FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV TZ America/Chicago

# Update the package list and install necessary packages
RUN apt update && \
    apt install software-properties-common -y && \
    apt update && \
    add-apt-repository ppa:deadsnakes/ppa -y && \
    apt update && \
    apt-get update && \
    apt-get install -y \
    build-essential \
    gcc \
    make \
    wget \
    curl \
    gzip \
    unzip \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    liblzma-dev \
    autotools-dev \
    autoconf \
    python3.12 \
    python3-pip \
    cython \
    git \
    bash
RUN apt-get install -y pkg-config

# install python environment with the super-fast uv solver
COPY requirements.txt requirements.txt
RUN python3 -m pip install uv
RUN python3 -m uv pip install -r requirements.txt

# configure LLM
RUN llm embed -m mini-l6 -c 'hello'
RUN llm sentence-transformers register \
    all-mpnet-base-v2 \
    --alias mpnet

# modify LLM source code in `llm_cluster.py` to have it accept floats
# RUN
