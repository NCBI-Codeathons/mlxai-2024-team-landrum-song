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
    python3.10 \
    python3.10-distutils \
    python3-pip \
    cython \
    git \
    bash
RUN apt-get install -y pkg-config python-html5lib

# directly install pip from pypa
RUN curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10

# make sure the image uses the correct version of Python
RUN cd /usr/bin && \
    unlink python && \
    ln -s /usr/bin/python3.10 python && \
    unlink python3 && \
    ln -s /usr/bin/python3.10 python3

# install python environment
RUN python3 -m pip install pandas polars biopython rich flytekit lxml llm

# configure LLM
RUN llm install llm-cluster
RUN llm install llm-sentence-transformers

# modify LLM source code in `llm_cluster.py` to have it accept floats
RUN rm /usr/local/lib/python3.10/dist-packages/llm_cluster.py
COPY assets/llm_cluster.py /usr/local/lib/python3.10/dist-packages/llm_cluster.py

# register the models
RUN llm sentence-transformers register all-mpnet-base-v2
RUN llm sentence-transformers register all-MiniLM-L12-v2
RUN llm sentence-transformers register multi-qa-mpnet-base-dot-v1

# install nextflow
RUN python3 -m pip install nextflow

# set cache and config directories for llm
ENV LLM_USER_PATH "/scratch/.llm/.config"
ENV HF_HOME "/scratch/.llm"

# Run a bash shell by default when the container starts
CMD ["/bin/bash"]
