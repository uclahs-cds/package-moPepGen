FROM blcdsdockerregistry/bl-base:1.0.0 AS builder

COPY . /opt/moPepGen

ARG PYTHON_VER=3.8.11
ARG BIOPYTHON_VER=1.79
ARG MATPLOTLIB_VER=3.3.1

RUN conda create -qy -p /usr/local\
    python==${PYTHON_VER}

RUN cd /opt/moPepGen/ && \
    pip install . biopython==${BIOPYTHON_VER} matplotlib==${MATPLOTLIB_VER}

# Deploy the target tools into a base image
FROM ubuntu:20.04
COPY --from=builder /usr/local /usr/local

ENV MPLCONFIGDIR=/opt/matplotlib
RUN mkdir $MPLCONFIGDIR && \
    chmod 777 $MPLCONFIGDIR

LABEL maintainer="Chenghao Zhu <ChenghaoZhu@mednet.ucla.edu>"
