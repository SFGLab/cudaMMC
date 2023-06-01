FROM nvidia/cuda:11.8.0-devel-ubuntu22.04 as base

WORKDIR /cudaMMC
RUN apt-get update && apt-get install -y --no-install-recommends \
    && apt-get install -y cmake ninja-build && rm -rf /var/lib/apt/lists/*
COPY . .
RUN mkdir build && cd build && cmake ../ -DCUDA_ARCH="70;75;80;86" -GNinja && ninja
