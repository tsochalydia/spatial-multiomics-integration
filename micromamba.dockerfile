# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:2.8

LABEL maintainer="Erkin <erkin.acar@bsse.ethz.ch>" \
      description="Micromamba 2.8 - spatial env, usable by any runtime UID/GID, called as an execution image by external Snakemake"

ARG PYTHON_VERSION=3.11

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    DEBIAN_FRONTEND=noninteractive

USER root

RUN apt-get update && apt-get install -y \
    build-essential gfortran pkg-config cmake \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    zlib1g-dev \
    libpng-dev libjpeg-dev libtiff5-dev \
    libcairo2-dev libfreetype6-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libblas-dev liblapack-dev libopenblas-dev libglpk-dev \
    libhdf5-dev \
    wget git libfftw3-dev libgsl-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# All three envs are built from pinned env files (exact package versions).
COPY --chmod=755 spatial_env.yml /tmp/spatial_env.yml
COPY --chmod=755 scvi_env.yml /tmp/scvi_env.yml
COPY --chmod=755 int_retrieval_env.yml /tmp/int_retrieval_env.yml

WORKDIR /tmp

RUN micromamba create -y -f /tmp/spatial_env.yml -n spatial \
    && micromamba clean --all --yes

RUN micromamba create -y -f /tmp/scvi_env.yml -n scvi_env \
    && micromamba clean --all --yes

RUN micromamba create -y -f /tmp/int_retrieval_env.yml -n int_retrieval_env \
    && micromamba clean --all --yes

USER root

# No env preloaded onto PATH - call the interpreter directly, or use
# `micromamba run -n spatial ...` if a rule needs full activation.
RUN chmod -R a+rX /opt/conda


COPY --chmod=755 entrypoint.sh /usr/local/bin/entrypoint.sh


ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]
