FROM ubuntu:22.04

ENV DEBIAN_FRONTEND noninteractive

# Install dependencies
RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y \
        python3 \
        clang \
        clang-tidy \
        g++ \
        git \
        gdb \
        make \
        cmake \
        cmake-curses-gui \
        mpi-default-dev \
        libblas-dev \
        liblapack-dev \
        libhdf5-mpi-dev \
        libboost-filesystem-dev \
        libboost-program-options-dev \
        libboost-test-dev \
        libboost-log-dev \
        bzip2 \
        sudo

RUN useradd -mr -u 10000 builduser
RUN adduser builduser sudo
RUN echo "builduser:builduser" | chpasswd

USER builduser
WORKDIR /home/builduser
