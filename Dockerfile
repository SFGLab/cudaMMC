FROM nvidia/cuda:11.8.0-devel-ubuntu22.04 as build
LABEL stage=build

WORKDIR /home/cudaMMC/cudaMMC
RUN apt-get update && apt-get install -y --no-install-recommends \
    && apt-get install -y cmake ninja-build && rm -rf /var/lib/apt/lists/*
COPY . .
RUN mkdir build && cd build && cmake ../ -DCUDA_ARCH="60;70;75;80;86" -GNinja && ninja

FROM nvidia/cuda:11.8.0-base-ubuntu22.04 as cudammc
WORKDIR /home/cudaMMC
RUN apt-get update && apt-get install -y --no-install-recommends vim python3 python3-pip && rm -rf /var/lib/apt/lists/*
RUN pip install numpy
COPY --from=build /home/cudaMMC/cudaMMC/build/cudaMMC .
COPY ./example_data ./example_data
COPY ./benchmark ./benchmark
