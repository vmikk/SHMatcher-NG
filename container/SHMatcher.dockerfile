
## Build stage 1 (Rust and Cargo)
FROM rust:1.89.0-slim AS rust
RUN cargo install runiq

## Build stage 2 - Build DuckDB + Exon extension
FROM rust:1.89.0-slim AS exonbuilder

RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
    build-essential cmake git curl ca-certificates pkg-config python3 \
    libssl-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

ENV RUSTUP_TOOLCHAIN=1.89.0-x86_64-unknown-linux-gnu
ENV CARGO=/usr/local/cargo/bin/cargo
ENV RUST_TOOLCHAIN=1.89.0-x86_64-unknown-linux-gnu
ENV PATH="/usr/local/cargo/bin:/usr/local/rustup/toolchains/${RUSTUP_TOOLCHAIN}/bin:${PATH}"

# Ensure the requested toolchain is installed and active
RUN rustup toolchain install 1.89.0-x86_64-unknown-linux-gnu \
  && rustup default 1.89.0-x86_64-unknown-linux-gnu \
  && rustup show \
  && which rustc \
  && rustc --version \
  && which cargo \
  && cargo --version

WORKDIR /opt
RUN git clone --depth 1 https://github.com/wheretrue/exon-duckdb.git
WORKDIR /opt/exon-duckdb
## Ensure submodules are up to date
RUN make pull || true
## Pin DuckDB submodule to match the official CLI version used at runtime
RUN cd duckdb \
  && git fetch --tags \
  && git checkout v1.3.2 \
  && cd ..
## Make DuckDB aware of this out-of-tree extension
RUN sh -lc 'printf "%s\n" "duckdb_extension_load(exon SOURCE_DIR /opt/exon-duckdb)" > duckdb/extension/extension_config_local.cmake'
## Build release (produces duckdb CLI and exon.duckdb_extension)
RUN make CLIENT_FLAGS="-DRust_TOOLCHAIN=${RUST_TOOLCHAIN} -DRust_COMPILER=/usr/local/rustup/toolchains/${RUSTUP_TOOLCHAIN}/bin/rustc -DRust_CARGO=/usr/local/cargo/bin/cargo" release

## Build stage 3 - Main
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
    ca-certificates locales procps \
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
    "minimap2>=2.30" \
  && ${CONDA_PREFIX}/bin/mamba clean --all --yes

## Add goclust, USEARCH, and DuckDB
RUN cd /opt/software \
    && wget https://github.com/vmikk/goclust/releases/download/0.3b/goclust \
    && chmod +x goclust \
    && mv goclust ${CONDA_PREFIX}/bin/ \
    && wget https://raw.githubusercontent.com/rcedgar/usearch_old_binaries/refs/heads/main/bin/usearch11.0.667_i86linux64 \
    && chmod +x usearch11.0.667_i86linux64 \
    && mv usearch11.0.667_i86linux64 ${CONDA_PREFIX}/bin/usearch \
    && wget https://github.com/duckdb/duckdb/releases/download/v1.3.2/duckdb_cli-linux-amd64.zip \
    && unzip duckdb_cli-linux-amd64.zip \
    && mv duckdb ${CONDA_PREFIX}/bin/duckdb \
    && rm duckdb_cli-linux-amd64.zip \
    && mkdir -p /usr/local/lib/duckdb/extensions

## Allow loading unsigned/local extensions (for exon)
ENV DUCKDB_ALLOW_UNSIGNED_EXTENSIONS=1

## Copy only the Exon extension (note that it is built against DuckDB v1.3.2)
COPY --from=exonbuilder /opt/exon-duckdb/build/release/extension/exon/exon.duckdb_extension /usr/local/lib/duckdb/extensions/exon.duckdb_extension

## Rust tools (from the Cargo-based stage)
COPY --from=rust /usr/local/cargo/bin/runiq /opt/software/conda/bin/runiq

WORKDIR /opt/software
