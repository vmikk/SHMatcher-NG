
## Build stage 1 (Rust and Cargo)
FROM rust:1.89.0-slim AS rust
RUN cargo install runiq

## Build stage 2 - Main
FROM debian:12-slim AS main

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV SHELL=/bin/bash
ENV CONDA_PREFIX="/opt/software/conda"
ENV PATH=${PATH}:"/opt/software/conda/bin/"

RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
    tar zip unzip pigz gzip zstd xz-utils bzip2 coreutils \
    curl wget git less gawk nano rename bc \
    ca-certificates locales  \
    libtre-dev libtre5 zlib1g zlib1g-dev liblzma-dev libbz2-dev libcurl4-openssl-dev libglpk-dev libglpk40 \
    build-essential \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

## Install conda
RUN mkdir -p /opt/software \
  && cd /opt/software \
  && curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
  && bash Miniforge3-Linux-x86_64.sh -u -b -p ${CONDA_PREFIX} \
  && rm Miniforge3-Linux-x86_64.sh \
  && ${CONDA_PREFIX}/bin/conda config --add channels bioconda \
  && ${CONDA_PREFIX}/bin/mamba update -y --all \
  && ${CONDA_PREFIX}/bin/mamba clean --all --yes

## Create conda environment and install software
RUN ${CONDA_PREFIX}/bin/mamba install -y \
    "vsearch>=2.30.0" \
    "seqkit>=2.10.1" \
    "bioawk" \
    "parallel>=20250622" \
    "csvtk>=0.31.0" \
    "ripgrep>=14.1.1" \
    "fd-find>=10.2.0" \
    "mmseqs2" \
  && ${CONDA_PREFIX}/bin/mamba clean --all --yes

## Add goclust and USEARCH
RUN cd /opt/software \
    && wget https://github.com/vmikk/goclust/releases/download/0.3b/goclust \
    && chmod +x goclust \
    && mv goclust ${CONDA_PREFIX}/bin/ \
    && wget https://raw.githubusercontent.com/rcedgar/usearch_old_binaries/refs/heads/main/bin/usearch11.0.667_i86linux64 \
    && chmod +x usearch11.0.667_i86linux64 \
    && mv usearch11.0.667_i86linux64 ${CONDA_PREFIX}/bin/usearch

## Rust tools (from the Cargo-based stage)
COPY --from=rust /usr/local/cargo/bin/runiq /opt/software/conda/bin/runiq

WORKDIR /opt/software
