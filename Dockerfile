FROM ghcr.io/seisscoped/container-base:ubuntu22.04_jupyterlab

COPY . .

RUN apt update && \
    apt install -y gcc gfortran make && \
    cd src && \
    make all && \
    cp ph2dt/ph2dt /usr/bin/ && \
    cp hypoDD/hypoDD /usr/bin/ && \
    cp hista2ddsta/hista2ddsta /usr/bin/

WORKDIR ${NB_HOME}
